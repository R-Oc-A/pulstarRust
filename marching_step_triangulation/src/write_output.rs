use crate::*;
use na::Vector3;

pub fn extract_points(triangles:&Tetrahedrization)->DataFrame{
    //extract point coordinates
    let points = &triangles.points;

    let coords:Vec<Vector3<f64>> = points.iter().map(|x|x.coords).collect();
    let point_ids:Vec<u32> = points.iter().map(|x|x.point_number as u32).collect();

    let coords_x:Vec<f64> = coords.iter().map(|x|x.x).collect();
    let coords_y:Vec<f64> = coords.iter().map(|x|x.y).collect();
    let coords_z:Vec<f64> = coords.iter().map(|x|x.z).collect();

    df![
        "point_id" => point_ids,
        "x coordinate" => coords_x,
        "y coordinate" => coords_y,
        "z coordinate" => coords_z,
    ].unwrap()
}


pub fn extract_triangles(triangle_collection:&Tetrahedrization)->DataFrame{
    let triangles = &triangle_collection.triangles;

    let first_vertex:Vec<u32> = triangles.iter().map(|x|
        x.vertices[0].point_number as u32
    ).collect();
    let second_vertex:Vec<u32> = triangles.iter().map(|x|
        x.vertices[1].point_number as u32
    ).collect();
    let third_vertex:Vec<u32> = triangles.iter().map(|x|
        x.vertices[2].point_number as u32
    ).collect();

    df![
        "first vertex" => first_vertex,
        "second vertex" => second_vertex,
        "third vertex" => third_vertex
    ].unwrap()
}

/*pub fn get_writer(file_name:&str)->ParquetWriter<>{
    let mut file=std::fs::File::create(file_name).expect("unable to create file");
    ParquetWriter::new(&mut File)    
}*/