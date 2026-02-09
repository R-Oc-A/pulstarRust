/*use tch::{kind::Kind, no_grad_guard, CModule, Device, IValue, Tensor};

fn regresor_templ() -> tch::Result<()> {
    let device = Device::cuda_if_available();
    let model = CModule::load_on_device("model_lines.ptrom", device)?;
    // The input is a 4 element tensor with Teff, logg, metallicity and vsini
    let input: Tensor = Tensor::rand(&[4], (Kind::Float, device));
    // We only care about inference so we can omit the gradients
    let _guard = no_grad_guard();
    // let output = model.forward_ts(&[input])?;
    // Calling model functions other than forward (see line above) is a less straightforward
    let output_iv = model.method_is("forward_ipca", &[IValue::from(input)])?;
    let output = match output_iv {
        IValue::Tensor(t) => t,
        other => panic!("forward_ipca returned non-tensor: {:?}", other),
    };
    println!("{:?}", output.size());
    let slice = output
        .slice(0, 0, 20, 1) // first 20 points just as an example
        .to_device(Device::Cpu) // at the moment this is redundant
        .to_kind(Kind::Float)
        .contiguous();
    let n = slice.numel();
    let mut vals = vec![0f32; n as usize];
    slice.copy_data(&mut vals, n);
    println!("{:?}", vals);
    Ok(())
}
*/

// grillas de nadya van de 3000 a 10000 angstroms

// 1 incluir pytorch to rust 
// 2 graficar los flujos
// 3 programar jalar el regresor al proyecto
// 4 construir los Array3 
// 5 generar las estructuras de datos [SpectralGrid] para las grillas de Pablo.


// Explorar Pyo3
