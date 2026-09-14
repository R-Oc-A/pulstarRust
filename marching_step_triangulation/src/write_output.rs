use crate::*;
use na::Vector3;


/// Module used to compute surface normal, centroid of a triangle..and other stuff.
pub mod triangle_functions;

#[derive(Clone)]
pub struct ExtractedTriangulation{
    pub triangles:Vec<[usize;3]>,
    pub points:Vec<Point>
}

pub fn extract_points(tetra:&Tetrahedrization)->DataFrame{
    //extract point coordinates
    let points = &tetra.points;

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


pub fn extract_triangles(tetra:&Tetrahedrization)->DataFrame{
    let triangles = &tetra.triangles;

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

pub fn extract_triangles_from_dataframe(points_df:&DataFrame,triangles_df:&DataFrame)->ExtractedTriangulation{

    // sort data points
    let points_sorted = points_df.clone().lazy().sort(["point_id"],Default::default()).collect().unwrap();


    //First with points
    let point_id = extract_column_as_vectorusize("point_id", &points_sorted);
    let coords_x = extract_column_as_vectorf64("x coordinate", &points_sorted);
    let coords_y = extract_column_as_vectorf64("y coordinate", &points_sorted);
    let coords_z = extract_column_as_vectorf64("z coordinate", &points_sorted);

    let mut point_collection:Vec<Point> = Vec::new();
    for (index,point_number) in point_id.into_iter().enumerate(){
        println!("extracting point collection");
        println!("point {}, index {}",point_number,index);
        let point = Point{
            coords:Vector3::from([
                coords_x[index],
                coords_y[index],
                coords_z[index]
            ]),
            point_number:point_number
        };
        point_collection.push(point);
    }

    //Second with triangles
    let first_vertices = extract_column_as_vectorusize("first vertex", triangles_df);
    let second_vertices = extract_column_as_vectorusize("second vertex", triangles_df);
    let third_vertices = extract_column_as_vectorusize("third vertex", triangles_df);
    let mut triangle_collection:Vec<[usize;3]> = Vec::new();
    for (index,first_vertex) in first_vertices.into_iter().enumerate(){
        let vertices = [
                first_vertex,
                second_vertices[index],
                third_vertices[index],
            ];
        triangle_collection.push(vertices);
    }

    ExtractedTriangulation { triangles: triangle_collection, points: point_collection }
    
}


/// This function takes a polars data frame and returns all of the values from a given column that holds f64 values. 
/// ### Arguments: 
/// * `column_name` - a string slice that holds the name of a column. The column should hold f64 values.
/// * `df`- a polars DataFrame
/// ### Returns:
/// * `Vec<f64>` - a vector that contains all of the values on the column.
fn extract_column_as_vectorf64(column_name: &str,df:&DataFrame)->Vec<f64>{
    let column = df.column(column_name).unwrap();
    let array = column.f64().unwrap();
    let extracted_vector:Vec<f64> = array.to_vec().into_iter().map(
        |x|match x{
            None => {panic!("error while extracting the data from a column")},
            Some(value )=>{value}
        }
    ).collect();
    extracted_vector
    //column.f64().unwrap().into_iter().flatten().collect()
}

/// This function takes a polars data frame and returns all of the values from a given column that holds f64 values. 
/// ### Arguments: 
/// * `column_name` - a string slice that holds the name of a column. The column should hold f64 values.
/// * `df`- a polars DataFrame
/// ### Returns:
/// * `Vec<usize>` - a vector that contains all of the values on the column.
fn extract_column_as_vectorusize(column_name: &str,df:&DataFrame)->Vec<usize>{
    let column = df.column(column_name).unwrap();
    let array = column.u32().unwrap();
    let extracted_vector:Vec<usize> = array.to_vec().into_iter().map(
        |x|match x{
            None => {panic!("error while extracting the data from a column")},
            Some(value )=>{value as usize}
        }
    ).collect();
    extracted_vector
    
}

/*pub fn get_writer(file_name:&str)->ParquetWriter<>{
    let mut file=std::fs::File::create(file_name).expect("unable to create file");
    ParquetWriter::new(&mut File)    
}*/

impl Tetrahedrization{
    pub fn triangulation_output(&mut self)->ExtractedTriangulation{
        let points_df = self.extract_points();
        let triangles_df = self.extract_triangles();

        extract_triangles_from_dataframe(&points_df, &triangles_df)
    }
}