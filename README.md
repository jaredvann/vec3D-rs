# vec3D

[![Build Status](https://travis-ci.org/jaredvann/vec3D-rs.svg?branch=master)](https://travis-ci.org/jaredvann/vec3D-rs)
[![Latest Version](https://img.shields.io/crates/v/vec3D.svg)](https://crates.io/crates/vec3D)
[![Documentation](https://docs.rs/vec3D/badge.svg)](https://docs.rs/vec3D)

A minimal 3D Vector library in Rust. Designed with a preference towards conventions from physics. Inspired by the [CLHEP](http://proj-clhep.web.cern.ch/proj-clhep/) Hep3Vector class.

# Conventions

This module uses the convention for describing spherical coordinates as used
in the physics community as follows:

 - **r** - radial distance
 - **theta** - polar angle
 - **phi** - azimuthal angle

And cylindrical coordinates:

 - **r** - radial distance
 - **phi** - angle
 - **z** - height along z-axis

 All angles are in radians.

## Examples

```rust
use vec3D::Vec3D;

fn main() {
    // Simple Initialisation
    let vec1 = Vec3D::new(1.0, 2.0, 3.0);
    println!("{}", vec1);               // Prints: "[1.0, 2.0, 3.0]"

    // Operator overloads for clean code
    let vec2 = Vec3D::new(3.0, 4.0, 5.0);
    let vec3 = vec1 + vec2;
    println!("{}", vec3);               // Prints: "[4.0, 6.0, 8.0]"

    let vec4 = vec3 * 2.0;
    println!("{}", vec4);               // Prints: "[8.0, 12.0, 16.0]"

    // Common vector operations
    let dot_product = vec3.dot(vec4);
    println!("{}", dot_product);        // Prints: "232"

    let vec5 = Vec3D::new(1.0, 0.0, 0.0);
    let vec6 = Vec3D::new(0.0, 1.0, 0.0);
    let cross_product = vec5.cross(vec6);
    println!("{}", cross_product);      // Prints: "[0.0, 0.0, 1.0]"

    // Plus initialisations from spherical/cylindrical coordinates, rotations and more
}
```

## License

This project is licensed under the MIT license ([LICENSE](LICENSE) or http://opensource.org/licenses/MIT)