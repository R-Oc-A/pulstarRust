use ndarray::{Array1,Array2,Axis};
use ndarray::prelude::*;
use ndarray_linalg::*;
use crate::type_def::PI;
use itertools::Itertools;


/// This function calculates Hough functions. It wont give them in the usual function that depends on cos(θ) but rather an equally spaced array of colatitude angles.
///
/// it is left to compute how to get the proper index. I think I will likely initialize the matrix that contains colatitude values and associated hough function values. 
/// This function also calculates the derivatives, which will be useful to compute the velocities. 
pub fn hough(
    q:f64,
    l:u16,
    m:i16,
    npts:usize,
    lmbd:f64,
    extra:bool,
)->(f64,//lambda eigenvalue
    Array1<f64>,// mu = cos(θ) values in the [-1,1] range
    Array1<f64>,// hough_radial(mu)
    Array1<f64>,// hough_latitudinal(mu)
    Array1<f64>,// hough_azimuthal(mu)
    Array1<f64>,// d hough_radial(mu)/d mu
    Array1<f64>,// d hough_latitudinal(mu)/d mu
    Array1<f64>,// d hough_azimuthal(mu)/d mu
    )
{

    //enforce an even number of points
    let m_size = npts/2;
    let npts = m_size*2;

    // define parity
    let parity = (l as i16-m)%2;

    // Calculate the interior/root points (mu_i = Cos((2i-1)Pi/2N)
    // Calculate the interior/root points (mu_i = cos(((2i-1)Pi)/2N)), where i = 1,...,N , N=total number of collocation points.
    //let n = na::Vector1::new(vec![0.0;m_size]);
    //Define the coefficients of the differential equation for the radial Hough function. 
    let mut mu_vec:Vec<f64>=Vec::with_capacity(m_size);
    let mut s_vec:Vec<f64>=Vec::with_capacity(m_size);
    let mut denom_vec:Vec<f64> = Vec::with_capacity(m_size);
    let mut coeffs2_vec:Vec<f64> = Vec::with_capacity(m_size);
    let mut coeffs1_vec:Vec<f64> = Vec::with_capacity(m_size);
    let mut coeffs0_vec:Vec<f64> = Vec::with_capacity(m_size);
    let q_sqrd = q.powi(2);
    for index in 0..m_size{
        mu_vec.push((PI/(npts as f64) * (index as f64 + 0.5)).cos());//by defining the cosine here this way, you avoid the singular points.
        s_vec.push( (1.0-mu_vec[index].powi(2)).sqrt() );
        denom_vec.push( 1.0 - q_sqrd * mu_vec[index].powi(2) );
        coeffs2_vec.push( s_vec[index].powi(2) );
        coeffs1_vec.push( -2.0 * mu_vec[index] * (1.0-q_sqrd)/(denom_vec[index].powi(2)) );
        coeffs0_vec.push(
            q*(m as f64)*(1.0+q_sqrd*mu_vec[index].powi(2))/denom_vec[index].powi(2)
            -(m.pow(2) as f64)/s_vec[index] );
    }

    let mut mu = Array1::from_vec(mu_vec);
    let mut s = Array1::from_vec(s_vec);
    let denom = Array1::from_vec(denom_vec);
    let coeffs2 = Array1::from_vec(coeffs2_vec);
    let coeffs1 = Array1::from_vec(coeffs1_vec);
    let coeffs0 = Array1::from_vec(coeffs0_vec);

    // define the parity factor
    let pf = m%2;

    let mut d2 = Array2::<f64>::zeros((m_size,m_size));
    let mut d1= Array2::<f64>::zeros((m_size,m_size));
    let mut d0 = Array2::<f64>::zeros((m_size,m_size));
    //This are for the other parity
    let mut d1_other= Array2::<f64>::zeros((m_size,m_size));
    let mut d0_other = Array2::<f64>::zeros((m_size,m_size));

    for i in 0..m_size{
        for j in 0..m_size{
            
            let j_index = (if extra {2 * j as i16 + parity}
                else{2 * j as i16 +1 - parity}) as f64;
            let cij =(PI*j_index/(npts as f64) * 
            ((npts+i) as f64 + 0.5 )).cos();
            let sij =(PI*j_index/(npts as f64) * 
            ((npts+i) as f64 + 0.5) ).sin();
            
            if pf.abs() == 1{
                d0[[i,j]]=cij*s[i];
                d1[[i,j]]=j_index * sij - cij*mu[i]/s[i];
                d2[[i,j]]=( -j_index.powi(2) * cij * s[i].powi(2)
                    - j_index * sij * mu[i] * s[i] - cij)/ s[i].powi(3);
            }else{
                d0[[i,j]]=cij;
                d1[[i,j]]=j_index * sij / s[i];
                d2[[i,j]]=j_index * (
                    mu[i] * sij - j_index * cij * s[i]
                ) / s[i].powi(3);
            }
            if extra {
                let cij = (PI* j_index/(npts as f64) * ((npts + i) as f64 + 0.5 )).cos();
                let sij = (PI* j_index/(npts as f64) * ((npts + i) as f64 + 0.5 )).sin();

                if pf.abs() == 1{
                    d0_other[[i,j]]=cij;
                    d1_other[[i,j]]=j_index * sij / s[i]
                }else{
                    d0_other[[i,j]]=cij * s[i];
                    d1_other[[i,j]]=j_index * sij - cij * mu[i] / s[i]
                };
            }
        }         
    }   

    d0 = d0.inv().unwrap();
    d1 = d1.dot(&d0);
    d2 = d2.dot(&d0);

    if extra{
        d0_other = d0_other.inv().unwrap();
        d1_other = d1_other.dot(&d0_other);
    }

    let full = 
        Array2::from_diag(&coeffs2).dot(&d2) +
        Array2::from_diag(&coeffs1).dot(&d1) +
        Array2::from_diag(&coeffs0);

    let (vals,vecs) = full.eig().unwrap();

    println!("Here are the eigenvalues homie");
    println!("eigs structure {:?}",vals.shape());
    println!("eigvecs structure {:?}",vecs.shape());

    // Remove complex eigenvalues and positive ones. 
    //I think this can be done faster using indexed_iter
    let mut whatev:Vec<usize> = Vec::with_capacity(m_size);
    for item in vals.indexed_iter(){
        if item.1.im() == 0.0 && item.1.re() < 0.0 {whatev.push(item.0)}
    }
    let sifted_vals = vals.select(Axis(0),&whatev);
    let sifted_vecs = vecs.select(Axis(1),&whatev);
    println!("Here are the real eigenvalues homie");
    println!("eigs structure {:?}",sifted_vals.shape());
    println!("eigvecs structure {:?}",sifted_vecs.shape());

    //Sort array 
    //in python is easier
    //look at this references https://lemonfold.io/posts/2023/rust/sorting_ndarray/
    //and https://github.com/rust-ndarray/ndarray/blob/master/examples/sort-axis.rs
    //Here bellow is the ordered index vals.
    let idx:Vec<usize> = sifted_vals.slice(s![..])//gets a reference to the data
                .to_owned().as_slice().unwrap()//Creates an owned copy of the slice and returns a slice view
                .into_iter()//creates an iterator of the data
                .enumerate()//gives an index for the iterator
                .sorted_by(|a,b|
                    a.1.re().partial_cmp(&b.1.re()).unwrap())//For each iterator orders the items by making a partial comparison. The fact that this works seems like magic to me because it smells like it's doing bubble sort.
                .map(|x|x.0).collect();//Maps the result of the comparison by given the permutation of the indices to get an ordered array.
    
    let mut eigenvals = sifted_vals.select(Axis(0),&idx);
    let mut eigenvecs = sifted_vecs.select(Axis(1),&idx);
    
    let ind_pos = eigenvals.slice(s![..])
            .iter()
            .enumerate()
            .fold(0usize, |x,b|{
                if (eigenvals[x].re() - lmbd).abs() < (b.1.re()-lmbd).abs(){x}
                else{b.0}
            });
    
    //Eigenvalue for the radial Hough function differential equation
    let eigenval = eigenvals[ind_pos].re();
    let mut hough_r = eigenvecs.column(ind_pos).map(|x|x.re());
    let norm = hough_r.norm();
    if (1.0-norm).abs()>std::f64::EPSILON{
        hough_r = hough_r.map(|x|x/norm);        
    }        

    //last point should be positive...
    if let Some(val) = hough_r.last_mut(){
        *val= val.abs()};
    

    // Compute latitudinal hough function
    let mut coeffs1_ht:Array1<f64> = Array1::zeros((m_size));
    let mut coeffs0_ht:Array1<f64> = Array1::zeros((m_size));
    for (index,denom_item) in denom.slice(s![..]).into_iter().enumerate(){
        coeffs1_ht[index] = - s[index].powi(2) / denom_item;
        coeffs0_ht[index] = - m as f64 * q * mu[index]/ denom_item;
    }
    let full_ht:Array2<f64> = Array2::from_diag(&coeffs1_ht).dot(&d1)+
        Array2::from_diag(&coeffs0_ht);
    
    let mut hough_t = full_ht.dot(&hough_r);
    for (index,ht) in hough_t.slice_mut(s![..]).into_iter().enumerate(){
        *ht = *ht/s[index];
    }

    // Compute azimuthal hough function
    let mut coeffs1_hp:Array1<f64> = Array1::zeros((m_size));
    let mut coeffs0_hp:Array1<f64> = Array1::zeros((m_size));
    for (index,denom_item) in denom.slice(s![..]).into_iter().enumerate(){
        coeffs1_hp[index] = q * mu[index] * s[index].powi(2) / denom_item;
        coeffs0_hp[index] = m as f64 * q * mu[index]/ denom_item;
    }
    let full_hp:Array2<f64> = Array2::from_diag(&coeffs1_hp).dot(&d1)+
        Array2::from_diag(&coeffs0_hp);

    let mut hough_p = full_hp.dot(&hough_r);
    for (index,hp) in hough_p.slice_mut(s![..]).into_iter().enumerate(){
        *hp = *hp/s[index];
    }

    //calculate the extra terms
    let mut hough_rp:Array1<f64> = Array1::zeros(m_size);
    let mut hough_tp:Array1<f64> = Array1::zeros(m_size);
    let mut hough_pp:Array1<f64> = Array1::zeros(m_size);

    if extra{
        hough_rp = d1.dot(&hough_r);


        let mut coeffs2_htp:Array1<f64>=Array1::zeros(m_size);
        let mut coeffs1_htp:Array1<f64>=Array1::zeros(m_size);
        let mut coeffs0_htp:Array1<f64>=Array1::zeros(m_size);

        for (index,denom_item) in denom.slice(s![..]).into_iter().enumerate(){
            coeffs2_htp[index] = - s[index].powi(4)/denom_item;
            coeffs1_htp[index] = s[index].powi(2) * (-m as f64 * q * mu[index] * denom_item 
                + mu[index] + mu[index].powi(3)*q.powi(2)
                - 2.0 * mu[index] * q.powi(2))/denom_item.powi(2);
            coeffs0_htp[index] = m as f64 * q 
                * ( 2.0 * mu[index].powi(4) * q.powi(2)
                - mu[index].powi(2) * q.powi(2) -1.0)/denom_item;
        }
        
        let full_htp = 
            Array2::from_diag(&coeffs2_htp).dot(&d2) +
            Array2::from_diag(&coeffs1_htp).dot(&d1) +
            Array2::from_diag(&coeffs0_htp);
        
        hough_tp = full_htp.dot(&hough_r);

        let mut coeffs2_hpp:Array1<f64>=Array1::zeros(m_size);
        let mut coeffs1_hpp:Array1<f64>=Array1::zeros(m_size);
        let mut coeffs0_hpp:Array1<f64>=Array1::zeros(m_size);

        for (index,denom_item) in denom.slice(s![..]).into_iter().enumerate(){
            coeffs2_hpp[index] = q * mu[index] * s[index].powi(4)/denom_item;
            coeffs1_hpp[index] = s[index].powi(2) * ( m as f64 * denom_item + q
                    + q.powi(3) * mu[index].powi(2)
                    - 2.0 * mu[index].powi(2) * q
                )/ denom_item.powi(2);
            coeffs0_hpp[index] = m as f64 * mu[index] * (
                    1.0 
                    - q.powi(2) * mu[index].powi(2) 
                    + 2.0 * q.powi(2) 
                    - 2.0 * q.powi(2) * mu[index].powi(2) )/ denom_item.powi(2);
        }
        
        let full_hpp = 
            Array2::from_diag(&coeffs2_hpp).dot(&d2) +
            Array2::from_diag(&coeffs1_hpp).dot(&d1) +
            Array2::from_diag(&coeffs0_hpp);
        
        hough_pp = full_hpp.dot(&hough_r);
    }

    // Append the symmetric terms.
    // This is necessary because the Laplace tidal differential operator 
    // is invariant under the transformation mu->-mu;
    // solving the equations only in a hemisphere.
    // You need to paste the symmetric parts of the solutions. 
    append_reflection(& mut mu, true);
    //append_reflection(& mut s, true);

    if parity.abs() == 1{
        append_reflection(& mut hough_r, true);
        append_reflection(& mut hough_t, false);       
        append_reflection(& mut hough_p, true);       
    }else{
        append_reflection(& mut hough_r,false); 
        append_reflection(& mut hough_t, true);       
        append_reflection(& mut hough_p,false);       
    }

    if extra{
        if parity.abs() == 1 {
            append_reflection(& mut hough_rp, true);
            append_reflection(& mut hough_tp, false);       
            append_reflection(& mut hough_pp, true);       
        }else{
            append_reflection(& mut hough_rp,false); 
            append_reflection(& mut hough_tp, true);       
            append_reflection(& mut hough_pp,false);       
        }
    }
    (
        eigenval,
        mu,
        hough_r,
        hough_t,
        hough_p,
        hough_rp,
        hough_tp,
        hough_pp
    )        
}

fn append_reflection(arr:& mut Array1<f64>,change_sign:bool){
    if change_sign{
        arr.append(Axis(0),
            arr.slice(s![..;-1])
            .to_owned()
            .mapv(|x| -x)
            .slice(s![..])).unwrap()
    }else{
        arr.append(Axis(0),
            arr.slice(s![..;-1])
            .to_owned()
            .slice(s![..])).unwrap()
    }
}