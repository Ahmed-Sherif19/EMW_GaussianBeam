# Gaussian Beam Propagation and Focusing Simulation

MATLAB implementation of Gaussian beam propagation using the Angular Spectrum Method, demonstrating wave optics and lens focusing.

## 🔬 Overview

Simulates Gaussian laser beam propagation through free space and focusing behavior with a thin lens. Developed for ECE242 Electromagnetic Waves course.

## ✨ Features

- **Angular Spectrum Method** for free space propagation using 2D Fourier transforms
- **Lens focusing simulation** with quadratic phase transformation
- **Comprehensive visualization** including:
  - 2D intensity maps
  - 1D cross-section at y=0
  - Phase distribution
  - 3D surface plots
- **Animation generation** showing complete beam evolution
- **High-quality exports** using `exportgraphics`


Check `outputs/` folder for generated figures and animation.

## 📊 Simulation Results

The simulation demonstrates beam evolution at key positions:

### Before Lens (Free Space Propagation)
- **z=0**: Perfect circular Gaussian shape at origin
- **z=z₀**: Beam expansion due to natural wave spreading  
- **z=2z₀**: Further spreading with lower center intensity

### After Lens Application
- **z=0.5f**: Beam begins to concentrate, converging behavior
- **z=f**: Maximum beam concentration at focal point
- **z=2f**: Beam diverges again after focus

## 🔧 Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Wavelength | 2 μm | Free-space wavelength |
| Beam Waist | 10 μm | Initial w₀ |
| Focal Length | 0.5z₀ | Lens focal length |
| Grid Size | 256×256 | Spatial sampling |
| Spatial Step | 1 μm | Grid spacing |

## 📁 Repository Structure

```
gaussian-beam-propagation/
├── src/
│   └── Gaussian_beam_simulation.m     # Main MATLAB script
├── docs/
│   └── Waves_Project_Report.pdf       # Technical report
├── outputs/
│   ├── figures/
│   │   ├── fig1_initial.png           # Initial beam
│   │   ├── fig2a_z0.png              # At z=z₀
│   │   ├── fig2b_2z0.png             # At z=2z₀
│   │   ├── fig3_0.5f.png             # After lens z=0.5f
│   │   ├── fig4_f.png                # At focal point
│   │   └── fig5_2f.png               # After focus z=2f
│   └── animation/
│       └── complete_beam_propagation.avi
└── README.md
```

## 🔬 Key Functions Used

- `fft2`, `ifft2` for Fourier transforms
- `meshgrid`, `imagesc`, `surf` for spatial grids and visualization
- `exportgraphics` for saving figures
- `VideoWriter` for generating animation

## 👥 Team

**Team 1 - Spring 2025**
- Ahmed Sherif Mohamed (23P0414)
- Ahmed Belal Abdelrahman (23P0007)  
- Amr Ahmed Zaki (23P0011)

**Supervisors:** Dr. Hussein Kotb, Eng. Mona Kamel

## 📚 Documentation

- [Technical Report](docs/Waves_Project_Report.pdf) - Complete analysis and results
- [Animation Demo](outputs/animation/complete_beam_propagation.avi) - Beam evolution visualization

## 📖 References

[1] D. M. Pozar, *Microwave Engineering*, 4th ed. Hoboken, NJ, USA: Wiley, 2011.  
[2] F. T. Ulaby and U. Ravaioli, *Fundamentals of Applied Electromagnetics*, 7th ed. Boston, MA, USA: Pearson, 2014.
