extern crate core;

use num::complex::Complex;
use plotters::prelude::*;

static W: f32 = 0.1;
static GAMMA : f32 = 0.01;

static D_12: f32 = 1.0;
static E_0: f32 = 1.0;
static E_1: f32 = 1.0;
static E_2: f32 = 1.1;
static C1_0: Complex<f32> = Complex::new(1.0, 0.0);

static HBAR: f32 = 1.0;
static MAX_TIME: f32 = 30.0;
static DELTA_TIME: f32 = 0.001;

static W_0: f32 = (E_2 - E_1) / HBAR;
static DELTA_W: f32 = W - W_0;
static OMEGA_R: f32 = D_12 * E_0 / (HBAR * 2.0);
static mut C2_0: Complex<f32> = Complex::new(0.0, 0.0);

struct GraphData {
    f1 : Vec<f32>,
    f2 : Vec<f32>
}


fn dc_1(t: f32, c2t: Complex<f32>, c1: Complex<f32>) -> Complex<f32> {
    -OMEGA_R * Complex::new(0.0f32, DELTA_W * t).exp() * c2t / Complex::new(0.0f32, 1.0f32)
         - Complex::new(0.0, GAMMA) * Complex::new(0.0, -E_1/HBAR * t) * c1
}

fn dc_2(t: f32, c1t: Complex<f32>, c2: Complex<f32>) -> Complex<f32> {
    -OMEGA_R * Complex::new(0.0f32, -DELTA_W * t).exp() * c1t / Complex::new(0.0f32, 1.0f32)
        -  Complex::new(0.0, GAMMA) * Complex::new(0.0, -E_2/HBAR * t) * c2
}

unsafe fn calculate() -> GraphData {
    let mut t = 0.0f32;
    let mut fc1 : Vec<f32> = Vec::new();
    let mut fc2 : Vec<f32> = Vec::new();
    let mut c1 = C1_0;
    let mut c2 = C2_0;
    fc1.push(c1.re);
    fc2.push(c2.re);
    while t <= MAX_TIME {
        t += DELTA_TIME;
        let c1p = c1;
        c1 += dc_1(t, c2, c1) * Complex::new(DELTA_TIME, 0.0f32);
        c2 += dc_2(t, c1p, c2) * Complex::new(DELTA_TIME, 0.0f32);
        fc1.push((c1 * c1.conj()).re);
        fc2.push( (c2 * c2.conj()).re);
    }
    return GraphData {
        f1: fc1,
        f2: fc2
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    unsafe { C2_0 = (1.0 - C1_0.powf(2.0f32)).sqrt(); }

    let res = unsafe { calculate() };
    plt(res, "plot.png")?;
    Ok(())
}

fn plt(data : GraphData, name: &str) -> Result<(), Box<dyn std::error::Error>>{
    let root = BitMapBackend::new(name, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(5u32)
        .x_label_area_size(30u32)
        .y_label_area_size(30u32)
        .build_cartesian_2d(0f32..MAX_TIME as f32, -0.1f32..1.1f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            (0..=(MAX_TIME/DELTA_TIME) as i32).map(|x| x as f32 * DELTA_TIME).map(|x| (x, data.f1[(x / DELTA_TIME) as usize])),
            &RED,
        ))?
        .label("C1")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .draw_series(LineSeries::new(
            (0..=(MAX_TIME/DELTA_TIME) as i32).map(|x| x as f32 * DELTA_TIME).map(|x| (x, data.f2[(x / DELTA_TIME) as usize])),
            &BLUE,
        ))?
        .label("C2")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    Ok(())
}

