/// Marching method for triangulation of surfaces using the algorithms presented in
/// "Hartmann, E. (1998). A marching method for the triangulation of surfaces. The Visual Computer, 14(3), 95-108."
/// 
/// The algorithm is described in 4 steps and here I follow the same logic; thus the main functions for the triangulations are named step0, ..., step4
use std::{f64::consts::PI, rc::Rc};
use nalgebra as na;
use na::{Unit,Vector3,Rotation3};
use polars::prelude::*;



/// Module used to extract data from a [Tetrahedrization]. 
mod write_output;


/// Module used to compute surface normal, centroid of a triangle..and other stuff.
mod triangle_functions;


/// Error threshold for the computation of the Newton rhapson method of landing a point close to the implicit surface into the surface.
const ERROR_THR:f64 = 1.0e-8;
/// Maximum number of steps after which the computation of the Newton rhapson method stops.
const NOT_CONVERGING:u16 = 1000;


/// A [FrontPoint] is a member of either a front polygon or a border polygon.
/// 
/// These points are the ones that will be covered by triangles, after which they are taken out of the structure and substitued by the external points of the triangulation that cover them. 
/// 
/// on each frontal point there's also a local coordinate system that has a surface normal vector and two tangent vectors orthonormal to it.
#[derive(PartialEq,Clone)]
pub struct FrontPoint{
    /// A [Rc<Point>] that contains the coordinates and point id of the front point. This data is shared with a [Triangle] that has this point as one of its vertex and with the whole [Tetrahedrization].
    point:Rc<Point>,
    /// Surface normal vector pointing outward of the implicit surface. 
    surface_normal:Unit<Vector3<f64>>,
    /// Normal vector tangent to the implicit surface and orthogonal to the surface normal. 
    t1:Unit<Vector3<f64>>,
    /// Normal tangent vector orthogonal to both the surface normal and t1. 
    t2:Unit<Vector3<f64>>,
    /// Front angles are used to determine how many triangles will be put around this point. To compute them you take into acount the two closest neighboring [FrontPoint]s
    front_angle:f64,
    /// Flag that states that the front angle should be recomputed when this point is new, or a neighbour of this point is new.
    angle_changed:bool,
    /// If this point is a member of a Border polygon rather than a front polygon, then it wont be covered by triangles. 
    border_point:bool,
}

/// A point on a surface. It is defined in x,y,z coordinates and has also a point identifier. These points are shared between [Triangle], [FrontPoint], and [Tetrahedrization].
#[derive(PartialEq,Clone)]
pub struct Point{
    /// * `coords` - the x,y,z coordinates of the point q close to the surface passed as a &[f64;3]
    coords:Vector3<f64>,
    /// * `point_number` - a [usize] that identifies the new point to be created.
    point_number:usize
}

/// This is the abstraction of triangles. It is defined as a collection of 3 vertex
/// where each vertex has a [Rc] to the coordinates [[f64];3]] and a [Rc] of the point number of the triangulation.
pub struct Triangle{
    vertices:[Rc<Point>;3],//point coordinates and point number
}

/// This is the collection of [Triangle]s that discretize the surface;
pub struct Tetrahedrization
    {
    /// A [f64] value that is approximately the length of an edge of a given triangle.  
    pub delta_t:f64,

    /// A [Vec] collection of [Triangle] that describe the tetrahedrization;
    pub triangles:Vec<Triangle>,
    
    /// The collection of [Point]s that are use to produce the mesh;
    pub points:Vec<Rc<Point>>,

    /// Front polygons are used to advance the triangulation. They are collection of points that get erased once they are covered with triangles. 
    /// The first [Vec<Point>] collection is called the _actual front polygon_ and it's the one that will be covered with triangles. Once it is finished
    /// the second front polygon becomes the actual front polygon and so on. 
    pub front_polygons:Vec<Vec<FrontPoint>>,
    
    /// How many points have already been created.
    pub current_number_points:usize,

    /// * `potential` - a closure f represented by |coords:&[f64;3]|->f64  that gives the definition of the surface implicitly as 
    /// f(x)=0'
    pub potential:Box<dyn Fn(&Vector3<f64>)->f64>,

    /// * `grad_potential` - a closure gf represented by |coords:&[f64;3]|->[f64;3] that gives the gradient of f; (and this should not be zero)
    pub grad_potential:Box<dyn Fn(&Vector3<f64>)->Vector3<f64>>
}

/// This function computes the "polar angle" refered in the 
fn compute_polar_angle(vec3: &Vector3<f64>)->f64{    
    let atan = (vec3.y/vec3.x).atan();
    if vec3.y>=0.0 {
        //first quadrant
        if vec3.x>=0.0{
            atan}
        // second quadrant
        else {
            PI + atan
        }
    }else{ 
        //third quadrant
        if vec3.x<0.0{
            PI + atan
        }
        //fourth quadrant
        else {
            2.0 * PI + atan
        }
    }
}

impl Tetrahedrization{
    /// Creates a new instance of a [Tetrahedrization].
    /// ### Arguments: 
    /// * `delta_t` - A [f64] value that is approximately the length of a side of a triangle.
    /// * `potential` - A closure that evaluates the function that implicitly defines the surface on a given (x,y,z) point.
    /// * `grad_potential` - A closure that evaluates the gradient of the potential. 
    /// ### Returns: 
    /// * This function returns a [Tetrahedrization] instance where most of the fields are empty. 
    pub fn new<F,Gf>(delta_t:f64,potential: F,grad_potential:Gf)->Self
    where 
        F:Fn(&Vector3<f64>)->f64 + 'static,
        Gf:Fn(&Vector3<f64>)->Vector3<f64> + 'static,
    {
        Tetrahedrization{
            delta_t:delta_t,
            triangles: Vec::new(), points: Vec::new(), front_polygons: Vec::new(), 
            current_number_points: 0usize, potential: Box::new(potential),
            grad_potential: Box::new(grad_potential)}
        
    }
}

impl Tetrahedrization{
    /// This function determines a surface [Point] p that is near a given point q in the vicinity of a surface. 
    /// q-p need not be exaclty perpendicular to the surface
    /// ### Arguments:
    /// ### Returns:
    /// * - a new instance of [Point]
    pub fn surface_point(
        & mut self,
        coords:Vector3<f64>,
        point_number:usize)-> FrontPoint
        {   

            let potential=&self.potential;
            let grad_potential = &self.grad_potential;
            let mut u_k = coords.clone();
            
            let mut error = 100.0;
            let mut counter = 0;
    
            let mut f_uk = potential(&u_k);
            if f_uk.abs() > 1.0e-8*self.delta_t
            {
            while error>ERROR_THR && counter<NOT_CONVERGING && f_uk.abs()>1.0e-10 {
                let f = potential(&u_k);
                let gf =grad_potential(&u_k);
                let norm_gf_pow2 = gf.norm().powi(2);
                if norm_gf_pow2 == 0.0{ 
                    println!("I break in a maximum");
                    break;}
                
                let u_k1:Vector3<f64>= &u_k + &gf *(-f/norm_gf_pow2);
                
                error = (&u_k - &u_k1).norm();
                
                counter +=1;
                f_uk = f;
                u_k = u_k1;
            }
            if counter>=NOT_CONVERGING{println!("f_uk is now {}",f_uk)}
            }else{println!("I'm here, what now? {}",f_uk)}


            let surface_normal:Unit<Vector3<f64>> = Unit::new_normalize(grad_potential(&u_k));
            
            let t1:Unit<Vector3<f64>> = Unit::new_normalize(
                if surface_normal.x>0.5 || surface_normal.y>0.5{
                    Vector3::from([surface_normal.y,-surface_normal.x,0.0])
                }else{Vector3::from([-surface_normal.z,0.0,surface_normal.x])});
            
            let t2:Unit<Vector3<f64>> = Unit::new_normalize(
                surface_normal.clone().cross(&t1));
            
            let new_point= 
                Point{
                    coords:u_k,
                    point_number:point_number
                };

            FrontPoint{
            point:Rc::new(new_point),
            surface_normal:surface_normal,
            front_angle:0.0,
            t1:t1,
            t2:t2,
            angle_changed:true,
            border_point:false,
        }
    }
}


impl Point{
    /// Computes the distance of Two points (Using the L2 norm.).
    fn distance(&self,point:&Point)->f64{
        (self.coords-point.coords).norm()
   }
   
}

impl Tetrahedrization{
    /// This function returns the length of the actual front polygon, that is, the number of points it contains. 
    fn afp_len(&self)->usize{
        self.front_polygons[0].len()
    }

    /// This function returns the indices of the left and right neighbours. 
    fn neighbors_indices(&self,point_index:usize)->(usize,usize){
        let afp_len = self.afp_len();
        ((point_index+afp_len -1)%afp_len,
        (point_index +1)% afp_len)
    }
}

//----------------------------------------
//-------------step 0---------------------
//----------------------------------------
/// This function begins the triangulation
/// ### Arguments:
/// * `potential` - a closure that implicitly defines the surface as f(x)=0. 
/// * `grad_potential` - a closure that is the gradient of f (and should not be 0)
/// * `starting_point` - a [[f64];3] that has the coordinates of a point close to the vicinity of a surface point. 
/// ### Returns:
/// * `FrontPolygon` - the first [Vec<Point>] collection that will be the actual front polygon
/// * `Triangulation` -the first 6 [Triangle]s of the triangulation
/// * `point_number` - the number id of the last point created at this stage.
pub fn step0(tetra: &mut Tetrahedrization,starting_point:Vector3<f64>){
    let front_point= tetra.surface_point(starting_point,0);

    let mut first_front_polygon:Vec<FrontPoint>=Vec::new();//vec![front_point.clone()];
    let mut first_six_triangles:Vec<Triangle>=Vec::new();
    let mut q:Vector3<f64>=Vector3::zeros();
    for i in 0..=5{
        let phase = i as f64 * PI/3.0;
        let cosp=phase.cos();
        let sinp = (1.0-cosp.powi(2)).sqrt()*phase.sin().signum();
        for j in 0..3usize{ 
            q[j] = front_point.point.coords[j] + tetra.delta_t*cosp*front_point.t1[j]
            + tetra.delta_t*sinp*front_point.t2[j];
        }
        let new_point = tetra.surface_point(q,i+1);
        first_front_polygon.push(new_point);
    }

    let mut first_points:Vec<Rc<Point>> = vec![Rc::clone(&front_point.point)];
    for j in 0..=5usize{
        let jp1= if j == 5 {0}else{j+1};
        let triangle = Triangle{
            vertices:[
            Rc::clone(&front_point.point),
            Rc::clone(&first_front_polygon[j].point),
            Rc::clone(&first_front_polygon[jp1].point)
            ]
        };
        first_six_triangles.push(triangle);
        first_points.push(
            Rc::clone(&first_front_polygon[j].point)
        );
    }

    tetra.front_polygons = vec![first_front_polygon];
    tetra.triangles =  first_six_triangles;
    tetra.points = first_points;
    tetra.current_number_points = 6usize;
}



//----------------------------------------
//-------------step 1---------------------
//----------------------------------------
impl Tetrahedrization{
///This function recalculates the actual front angle of a new point (or its neighbors) of the actual front polygon.
/// It must receive a new [Point] or
/// a neighbour of a new [Point]. The specification of left and right neibours must be done by the calling method or function. 
/// ### Arguments:
/// * `left_neighbour`- A reference to P_0i-1 if i>1 or P_0N if i=0
/// * `right_neighbour`- A reference to P_0i+1 if i<N or P_01 if i=N
/// ### Returns:
/// * This method returns nothing. It just updates the actual front angle.
    pub fn step1(&mut self, point_index: usize){
        let front_point = self.front_polygons[0][point_index].clone();
        if front_point.angle_changed{
            let (v1_index,v2_index) = self.neighbors_indices(point_index);
            
            let v1= &self.front_polygons[0][v1_index];
            let v2= &self.front_polygons[0][v2_index];

            let front_angle = Self::get_angle_between_points(&front_point, v1, v2);
            self.front_polygons[0][point_index].front_angle = front_angle;
            self.front_polygons[0][point_index].angle_changed = false;
        }
    }

    fn get_angle_between_points(front_point:&FrontPoint,v1:&FrontPoint,v2:&FrontPoint)->f64{

        let w1 = compute_polar_angle(&v1.project_in_tangent_space(front_point));       
        let w2 = compute_polar_angle(&v2.project_in_tangent_space(front_point));
        
        if w2>=w1{
            let w = w2-w1;
            w
        }else{
            let w=w2-w1+2.0*PI;
            w
            
        }
    }
}



//----------------------------------------
//-------------step 2---------------------
//----------------------------------------
impl Tetrahedrization {
    
    // Remarks 1, 2, and 4.
    /// Checks the distances of each point of the actual front polygon and either splits it or joins it with another front polygon. 
    /// ### Arguments:
    /// * `points_to_ignore` - a [usize] collection of point's id that must be left out of this iteration of a distance check. 
    /// * `point_index` - a [usize] indicating the current point of the actual front polygon where we check distances. 
    /// ### Returns:
    /// * `point_number` - a [usize] indicating the number of the point to be left out of further distance checks. 
    /// * `restart_distance_check` - a [bool] that indicates a change in the actual front polygon, so that distance checks must be done again.
    pub fn step2(& mut self,points_to_ignore:&[usize],
        point_index:usize)->((usize,usize),bool){
        //check if current point must be ignored.
        if !self.not_in_points_to_ignore(point_index, points_to_ignore){return ((0,0),false)}
        //select actual front polygon
        let afp = &self.front_polygons[0];
        let afp_len = self.afp_len();
        if afp_len <5{return ((0,0),false);}
        //check local neighbours of the actual front polygon
        let (v1_index,v2_index) = self.neighbors_indices(point_index);
        let (vv1_index,_) = self.neighbors_indices(v1_index);//v1_index -1;
        let (_,vv2_index) =self.neighbors_indices(v2_index);//v2_index+1;
        let mut break_i = false;
        let mut p0s:(usize,usize) = (0,0);
        let mut index_j =0usize;
        for i in point_index..afp_len{
            if i != point_index &&
            i != v1_index &&
            i != v2_index &&
            i != vv1_index &&
            i != vv2_index 
            {
                if afp[point_index].point.distance(
                    &afp[i].point)< self.delta_t{
                        break_i = true;
                        index_j=i;
                        p0s=(afp[point_index].point.point_number,
                            afp[i].point.point_number);
                        break;
                    }
            };
        }
            
        if break_i {
            self.front_polygons[0][point_index].angle_changed = true;
            self.front_polygons[0][index_j].angle_changed = true;
            self.split_actual_front_polygon(point_index,index_j);
            return (p0s,true)
        }
            
        //check for among other actual front polygons
        if self.front_polygons.len() !=1{
            let mut index_m=0usize;
            let mut break_m =false;
            for (m,front_polygon) in self.front_polygons.iter().enumerate(){
                if break_m {break;}
                index_m = m;
                if index_m != 0{
                    for (j,front_point) in front_polygon.iter().enumerate(){
                        index_j=j;
                        if points_to_ignore.contains(&front_point.point.point_number)
                        {
                            continue};
                        
                        if afp[point_index].point.distance(
                            &front_point.point)<self.delta_t{
                            break_m=true;
                            p0s = (afp[point_index].point.point_number,
                                front_point.point.point_number);
                            break;
                        }
                    }
                }
            }
        
            if break_m{
                if !(self.checkbadpoints(
                    point_index,index_m,index_j)) //if not a bad point
                {   
                    self.front_polygons[0][point_index].angle_changed=true;
                    self.front_polygons[index_m][index_j].angle_changed=true;
                    self.join_front_polygons(
                    index_m, 
                    point_index, 
                    index_j);
                    return(p0s,true)
                }else{
                println!("I'm here, found some bad points, I guess.");
                return (p0s,true)}
            }
        }
        (p0s,false)
    } 

    
    /// This function  Joins the actual front polygon to another front polygon 'm' from the collection
    /// at some points 0i in the actual front polygon and mj in the other front polygon. 
    /// the new actual front polygon goes like this
    /// [p00, .., p0i, pmj, .., pmN, pm0, .., pmj, p0i, .., p0N]
    /// the points p0i and pmj appear twice but after some operations the first time they appear is deleted. 
    fn join_front_polygons(&mut self,
        polygon_index:usize,
        index_0i:usize,
        index_mj:usize)
        {
        let extracted_polygon=self.front_polygons.remove(polygon_index);

        //make new actual front polygon
        let actual_front_polygon = self.front_polygons[0].clone();
        
        let (first_slice,_) = actual_front_polygon.split_at(index_0i+1);
        let (_,second_slice) = extracted_polygon.split_at(index_mj);
        let (third_slice,_) = extracted_polygon.split_at(index_mj+1);
        let (_,fourth_slice)=actual_front_polygon.split_at(index_0i);
        let new_actual_front_polygon = 
        Vec::from(
            [
                [first_slice,second_slice].concat(),
                [third_slice,fourth_slice].concat()
            ].concat()
        );
        self.front_polygons[0]=new_actual_front_polygon;

        //recalculate frontal angles
        self.front_polygons[0][index_0i].angle_changed=true;
        self.front_polygons[0][index_0i+1].angle_changed=true;
        
        self.step1(index_0i);
        self.step1(index_0i+1);



        let p0i = self.front_polygons[0][index_0i].clone();
        let p0j = self.front_polygons[0][index_0i+1].clone();

        //Triangulate p0i and p0j beginning with the one that has the smallest front angle
        if p0i.front_angle<=p0j.front_angle{
            self.step3(index_0i);

            let mut index_0j=index_0i;
            for i in 0..self.afp_len(){
                if self.front_polygons[0][i].point.point_number==p0j.point.point_number{
                    index_0j = i;
                    break
                }
            }
            self.step1(index_0j);
            self.step3(index_0j);
             
            
        }else{
            self.step3(index_0i+1);
            
            // update the front angle
            self.step1(index_0i);
            self.step3(index_0i);
        }
    }

    /// This function splits the actual front polygon at a given index. The new front polygon
    /// The new front polygon goes from i to j
    /// and the new actual front polygon is a stitch from 0 to i and from j to N;
    fn split_actual_front_polygon(&mut self, index_i:usize, index_j:usize){
    // Make new front polygon:
        let mut new_actual_front_polygon = self.front_polygons[0].clone();
        let new_front_polygon:Vec<FrontPoint> = new_actual_front_polygon.drain(index_i..=index_j).collect();
        // Make new actual front polygon.
        let p0i = new_front_polygon[0].clone();
        let p0j = new_front_polygon.last().unwrap().clone();
        new_actual_front_polygon.insert(index_i,p0i.clone());
        new_actual_front_polygon.insert(index_i+1,p0j.clone());
    
        self.front_polygons[0]=new_actual_front_polygon;
        self.front_polygons.push(new_front_polygon);
    }

    //Check for bad points.
    fn checkbadpoints(& self, index_0i:usize,index_m:usize,index_j:usize)->bool{
        let (v1_index,_) = self.neighbors_indices(index_0i);
        let front_point = &self.front_polygons[0][index_0i];
        let v1 = &self.front_polygons[0][v1_index];
        let v2 = &self.front_polygons[index_m][index_j];
        Self::get_angle_between_points(front_point, v1, v2) >= front_point.front_angle
    }
}

impl FrontPoint{
    /// This function gives the coordinates of a [Point] in the [tangent space](https://en.wikipedia.org/wiki/Tangent_space) of another [Point]
    /// ### Arguments:
    /// * `reference_point`- a reference to a [Point] that contains the local orthonormal system defined by two tangent vectors and a normal.
    /// ### Returns:
    /// * `new_coords` - a [[f64];3] collection of numbers that have the x',y',z' coords.
    fn project_in_tangent_space(&self,reference_point: &FrontPoint)->Vector3<f64>{
        // move vector into reference point self.coords-point.coords

        let moved = self.point.coords-reference_point.point.coords;
        // do the dot product with t1,t2,normal, which are the new coords\
        let xi = moved.dot(&reference_point.t1.into_inner());
        let eta = moved.dot(&reference_point.t2.into_inner());
        let zeta = moved.dot(&reference_point.surface_normal.into_inner());

        // return this coords in the tangent space
        Vector3::<f64>::from([xi,eta,zeta])
    }

        }



impl Tetrahedrization{

    /// Complete the triangulation at point of the actual front polygon.
    /// ### Arguments:
    /// * `point_index` - a [usize] that states which [FrontPoint] of the actual front polygon is to be triangulated. 
    /// ### Returns:
    /// This function does not return anything 
    fn step3(&mut self,
        point_index:usize,)
        {
        // Determine the neighbors
        let (v1_indx,v2_indx) = self.neighbors_indices(point_index);
        let v1= self.front_polygons[0][v1_indx].clone();
        let v2 = self.front_polygons[0][v2_indx].clone();
        let front_point = self.front_polygons[0][point_index].clone();
        // Determine the number of triangles nt to be generated and the interior angles dw of each of the new triangles.
        let (nt,dw)=step3_2_number_of_triangles(self.delta_t,&front_point,&v1,&v2);
        // get new triangles
        let new_points = 
        self.generate_triangles(
            nt, dw,
            &front_point, &v1, &v2,
            );
        // Update actual_front_polygon
        self.step3_4_renew_actual_front_polygon(
            point_index, 
            new_points);
    }
}


///This function determines the number of triangles nt to be generated 
/// ### Arguments:
/// * point - a reference to a [RefCell<Point>] that will be used for the new triangles
/// * v1,v2 - a reference to the [RefCell<Point>] neighbors of point.
/// ### Returns:
/// (nt,dw)- a tuple of two [f64] values where the first one is nt and the second one is the angle of the triangle that has the reference point as its vertex.
fn step3_2_number_of_triangles(delta_t:f64,front_point:&FrontPoint,v1:&FrontPoint,v2:&FrontPoint)->(u16,f64){
    let w = front_point.front_angle;
    let mut nt = (3.0*w/PI).floor() as u16 +1;
    let mut dw = w/(nt as f64);
    // correct dw for extreme cases
    if dw<=0.65 && nt>1{
        nt-=1;
        dw = dw/(nt as f64);
    }else if nt==1 && dw>0.65 && v1.point.distance(&v2.point)>1.2 * delta_t{
        nt = 2;
        dw = 0.5*w;
    }else if w<3.0 && (
        v1.point.distance(&front_point.point)<=0.5 * delta_t ||
        v2.point.distance(&front_point.point)<=0.5 * delta_t){
            nt=1;
            dw = w;
    }
    (nt,dw)
}


impl Tetrahedrization{

/// This function generates the new triangles that surround a given [FrontPoint]
/// ### Arguments:
/// * `nt` - A [u16] that indicates the number of triangles that will sorrownd a given [FrontPoint]
/// * `dw` - A [f64] value that indicates the internal angle of the new triangles. 
/// * `front_point` - The [FrontPoint] that will be covered by triangles
/// * `v1` - The first neighbor of front_point
/// * `v2` - The second neighbor of the front_point. 
/// ### Returns:
/// This function returns a [Option] with the following variants:
/// * [Some<Vec<FrontPoint>>] - In case new [FrontPoint]s where inserted to create the new triangles. This points are going to be used to update the actual front polygon. 
/// * [None] - In case there's no new points.
fn generate_triangles(
    & mut self,
    nt:u16,dw:f64,
    front_point:&FrontPoint,
    v1:&FrontPoint,
    v2:&FrontPoint)
    ->Option<Vec<FrontPoint>>
    {
    if nt==1{
        let triangle = Triangle{
            vertices:[
                Rc::clone(&front_point.point),
                Rc::clone(&v1.point),
                Rc::clone(&v2.point),
            ]
        };
        self.triangles.push(triangle);
        None
    }else{
    
    let mut q_0:Vector3<f64> = v1.project_in_tangent_space(front_point);
    q_0.z =0.0;
    q_0 = Unit::new_normalize(q_0).into_inner();
    q_0 *= self.delta_t;
    //let q_0:Vector3<f64> = (Unit::new_normalize(v1.project_in_tangent_space(front_point))).into_inner() * self.delta_t;

    let mut q_i:Vec<Vector3<f64>>=Vec::with_capacity(nt as usize-1);
    //define the new points in the tangent space
    for i in 1..nt{
        let angle = i as f64 *dw;
        //make rotation matrix
        let axisangle = Vector3::z() * angle;//rotate around z which is the surface normal by an angle
        let rot = Rotation3::new(axisangle);
        //make q_i to be added
        q_i.push(
            rot * q_0
        );
    }
    //convert the tangent space points into coordinates in the regular reference frame.
    q_i = q_i.iter().map(|point_qi|->Vector3<f64>{
        &front_point.point.coords +
        (&front_point.t1.into_inner() * point_qi.x)+
        &front_point.t2.into_inner() * point_qi.y}).collect(); 
        
    //get the new points 
    let mut p_i:Vec<FrontPoint> = vec![v1.clone()];
    for q in q_i.into_iter(){
        self.current_number_points +=1;
        p_i.push(self.surface_point(q,self.current_number_points));
    }
    p_i.push(v2.clone());
    
    for p in p_i.iter_mut(){
        p.angle_changed=true;
    }
    // define the new triangles
    for index in 0usize..p_i.len()-1{
        p_i[index].angle_changed=true;
        p_i[index+1].angle_changed=true;
        self.triangles.push(
            Triangle { vertices: [
                Rc::clone(&p_i[index].point),
                Rc::clone(&p_i[index+1].point),
                Rc::clone(&front_point.point)
                ] }
        )
    }
    // Return new points
    Some(p_i.drain(1..p_i.len()-1).collect())
    }
}

/// This function updates the actual front polygon
/// ### Arguments:
/// * `index_of_point_to_remove` - A [usize] that selects a point from the actual front polygon that will be removed. Some [FrontPoint]s are repeated in the array so it is important to find them using these indices.
/// * `new_points_to_insert` - An [Option] that contains the new points created while generating the new [Triangle]s.
fn step3_4_renew_actual_front_polygon(&mut self,
    index_of_point_to_remove:usize,
    new_points_to_insert:Option<Vec<FrontPoint>>){
    
    // update neighbours:
    
    let (v1_indx,v2_indx) = self.neighbors_indices(index_of_point_to_remove);
    self.front_polygons[0][v1_indx].angle_changed=true;
    self.front_polygons[0][v2_indx].angle_changed=true;
    if let Some(new_front_points)=new_points_to_insert{
        // Add new points to the point collection
        let mut new_simple_points:Vec<Rc<Point>> = new_front_points.clone().into_iter().map(
            |front_point|->Rc<Point>{
                Rc::clone(&front_point.point)
            }
        ).collect();

        self.points.append(&mut new_simple_points);

        // update the actual front polygon
        let(first_slice,second_slice) = self.front_polygons[0].split_at(index_of_point_to_remove);
        let first_half = [first_slice,&new_front_points[..]].concat();
        let complete = [&first_half[..],&second_slice[1..]].concat();
        self.front_polygons[0]=complete;

    }else{
        self.front_polygons[0].remove(index_of_point_to_remove);
    }

}

}

//aa

//----------------------------------------
//-------------step 4---------------------
//----------------------------------------
impl Tetrahedrization{
    /// Following the marching algorithm reference, this function iterates steps 1 through 3 until the actual
    /// front polygon consist of only three points that generate a new triangle.
    /// If there is another (nonempty) front polygon left, it becomes the new actual front polygon and steps 1 through 3 are repeated. 
    /// If there are no more front polygons, then the triangulation is finished. 
    pub fn step4(&mut self){
        let mut counter = 0u32;
        while self.front_polygons.len()>0{
        //while counter < 4{
            counter += 1;
            // Here we add triangles and points to the actual front polygon until we're left with 3 points.
            if self.front_polygons[0].len()==3{
                println!("putting last triangle");
                self.put_last_triangle();
                println!("front polygons to be processed are {}",self.front_polygons.len());
                if self.front_polygons.len()>0{
                    println!("length of actual front polygon is {}",self.afp_len());
                    println!("--------------------");
                    continue;
                }else{
                    println!("----------------------------------------");
                    println!("---------Triangulation Finished---------");
                    println!("----------------------------------------");
                    continue;}
            }else if self.front_polygons[0].len()<3{
                println!("The length of the actual front polygon is {} so I don't know how to resume.",self.afp_len());
                break
            }else{
                if counter > 70 {
                    println!("69 iterations and it doesn't stop");
                    break
                }
                self.process_actual_front_polygon();
            }
        }
    }

    fn find_minimal_front_angle(&mut self)->usize{
            let (index_0m,_angle)=self.front_polygons[0]
            .iter()
            .enumerate()
            .fold((0usize,2.0*PI),|(indx_0,acc),(indx_0m,point)|{
                if point.front_angle<=acc{(indx_0m,point.front_angle)}
                else{(indx_0,acc)}
            });
        index_0m
    }

    fn not_in_points_to_ignore(& self,index:usize,points_to_ignore:&[usize])->bool{
        !points_to_ignore.contains(&self.front_polygons[0][index].point.point_number)
    }

    fn put_last_triangle(&mut self){
        if self.afp_len()==3{
            let afp = self.front_polygons.swap_remove(0);
            let triangle = Triangle{
            vertices:[
                Rc::clone(&afp[0].point),
                Rc::clone(&afp[1].point),
                Rc::clone(&afp[2].point),
            ]};
            self.triangles.push(triangle);
        }else{println!("why do you want to finish it?")}
    }

    fn process_actual_front_polygon(&mut self){
        let mut points_to_ignore = vec![0usize];
        //for _i in 0..2000{
        //while self.afp_len()>3 && self.afp_len()<200 {
        while self.afp_len()>3{
            for index_i in 0..self.afp_len(){
                self.step1(index_i);
            }
            let index_0m= self.find_minimal_front_angle();
            if self.front_polygons[0][index_0m].front_angle < 1.0
            {
                self.step3(index_0m);
                if self.afp_len()==3{
                    println!("breaking while finishing points with small front angles");
                    break
                }else if self.afp_len()<3{
                    println!("breaking badly while finishing points with small front angles");
                    println!("afp len is {}",self.afp_len());
                    break;
                }
                else{continue;}
            }
            let mut afp_changed =false;
            
            let mut new_points_to_ignore: (usize,usize) = (0,0);
            for index_i in 0..self.afp_len(){
                (new_points_to_ignore,afp_changed)=
                    self.step2(&points_to_ignore,index_i);        
                if afp_changed ||
                (new_points_to_ignore.0!=0&& new_points_to_ignore.1!=0)
                 {
                    points_to_ignore.push(new_points_to_ignore.0);
                    points_to_ignore.push(new_points_to_ignore.1);
                    break
                }
            }
            if afp_changed{
                if self.afp_len()<=3{
                    println!("breaking here after step2");
                    break;
                }
                else{continue}
            }
            if self.afp_len() == 3{
                println!("breaking here after distance checks");
                break
            }
            
            let index_0m = self.find_minimal_front_angle();
            self.step3(index_0m);
            
        }
    }

}


impl Tetrahedrization{
    /// Extracts the collection of [Point]s from a [Tetrahedrization] and returns it as a [DataFrame].
    /// ### Returns:
    /// * A Dataframe with the following columns
    /// * `|Point id|x coordinate|y coordinate|z coordinate|` - where the (x,y,z) are [f64] values and point id is a [u32].
    pub fn extract_points(&mut self)->DataFrame{
        write_output::extract_points(self)
    }

    /// Extracts the collection of [Triangle]s from a [Tetrahedrization] and returns it as a [DataFrame].
    /// ### Returns:
    /// * A Dataframe with the following columns
    /// * `|first vertex|second vertex|third vertex|` - where the vertices are [u32] values that correspond to a point id of a [Point].
    pub fn extract_triangles(&mut self)->DataFrame{
        write_output::extract_triangles(self)
    }
}