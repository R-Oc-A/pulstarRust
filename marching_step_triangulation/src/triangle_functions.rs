use super::*;

impl Triangle{
    /// By definition, the centroid of a triangle is located at the intersection of 2 medians; Thus it is the point located at 2/3 of the median between a point and a line. 
    /// ### Returns:
    /// * This function returns a [Point] with the coordinates of the location of the centroid of the triangle. 
    pub fn centroid(&self)->super::na::Vector3<f64>{
        //reference triangle
        //     c
        //   / | \
        //  /  *  \
        // a---d---b

        // First let's zet the center of coordinates around a
        let first_vertex = self.vertices[0].coords;//point a
        let second_vertex = self.vertices[1].coords;//point b
        let third_vertex =self.vertices[2].coords;// point c

        // Locate the mid point of a segment
        let ab = second_vertex - first_vertex;
        let d = 0.5*ab;

        // define the segment cd
        let c = third_vertex-first_vertex;
        let cd = d-c;

        //define the centroid
        let centroid = 0.3333333333333330 *cd; 

        //return the centroid in the original reference system
        centroid + first_vertex
    }

    /// This function computes the area and surface normal of a [Triangle] using the vector product
    pub fn area_and_surface_normal(&self)->(super::na::Vector3<f64>,f64){
        //reference triangle
        //    c
        //   / \
        //  a---b

        // First let's zet the center of coordinates around a
        let first_vertex = self.vertices[0].coords;//point a
        let second_vertex = self.vertices[1].coords;//point b
        let third_vertex =self.vertices[2].coords;// point c

        // define the segments
        let ab = second_vertex - first_vertex;
        let ac = third_vertex - first_vertex;

        // compute the area
        let cross_product = ab.cross(&ac);
        let area = cross_product.norm();
        let surface_normal = cross_product/area;
        
        //return area and surface normal
        (surface_normal,area*0.5)
    }

}

//to do unit testing.