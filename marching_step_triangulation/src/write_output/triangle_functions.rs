use super::*;
use std::ops::{Add, AddAssign, Mul};


impl ExtractedTriangulation{
    /// By definition, the centroid of a triangle is located at the intersection of 2 medians; Thus it is the point located at 2/3 of the median between a point and a line. 
    /// ### Returns:
    /// * a new [Point] with the coordinates of the location of the centroid of the triangle. 
    pub fn centroid(&self,index:usize)->super::na::Vector3<f64>{
        //reference triangle
        //     c
        //   / | \
        //  /  *  \
        // a-------b

        let first_vertex = self.points[self.triangles[index][0]].coords.clone();//point a
        let second_vertex = self.points[self.triangles[index][1]].coords.clone();//point b
        let third_vertex = self.points[self.triangles[index][2]].coords.clone();//point c

        //define the centroid
        let centroid = 0.3333333333333330 *first_vertex+second_vertex+third_vertex; 

        //return the centroid in the original reference system
        centroid
    }

    /// This function computes the area and surface normal of a [Triangle] using the vector product
    pub fn area_and_surface_normal(&self,index:usize)->(super::na::Vector3<f64>,f64){
        //reference triangle
        //    c
        //   / \
        //  a---b

        let first_vertex = self.points[self.triangles[index][0]].coords.clone();//point a
        let second_vertex = self.points[self.triangles[index][1]].coords.clone();//point b
        let third_vertex = self.points[self.triangles[index][2]].coords.clone();//point c

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

    /// This function generates a linear interpolator for a given triangle on it's centroid.
    pub fn interpolate_in_centroid<T>(&self,value_at_vertices:&[T;3])->T
    where
        T: Add<T,Output=T> + Mul<f64,Output = T> + AddAssign + Clone,
        f64: Mul<T,Output = T> 
    {   
        let mut interpolation = value_at_vertices[0].clone();

        for i in 1..=2{
            interpolation += value_at_vertices[i].clone()
        }
        interpolation*0.3333333333333330
    }
}

//TO DO unit testing.