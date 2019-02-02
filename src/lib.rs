// Copyright 2019 Jared Vann

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

//! A minimal 3D Vector library in Rust. Designed with a preference towards
//! conventions from physics. Inspired by the CLHEP Hep3Vector class.
//!
//! # Conventions
//!
//! This module uses the convention for describing spherical coordinates as used
//! in the physics community as follows:
//!
//!  - **r** - radial distance
//!  - **theta** - polar angle
//!  - **phi** - azimuthal angle
//!
//! ![spherical-coordinates-diagram](https://upload.wikimedia.org/wikipedia/commons/4/4f/3D_Spherical.svg)
//!
//! And cylindrical coordinates:
//!
//!  - **r** - radial distance
//!  - **phi** - angle
//!  - **z** - height along z-axis
//!
//!  All angles are in radians.
//!
//! # Examples
//!
//! ```
//! use vec3D::Vec3D;
//!
//! fn main() {
//!     // Simple Initialisation
//!     let vec1 = Vec3D::new(1.0, 2.0, 3.0);
//!     println!("{}", vec1); // Prints: "[1.0, 2.0, 3.0]"
//!
//!     // Operator overloads for clean code
//!     let vec2 = Vec3D::new(3.0, 4.0, 5.0);
//!     let vec3 = vec1 + vec2;
//!     println!("{}", vec3); // Prints: "[4.0, 6.0, 8.0]"
//!
//!     let vec4 = vec3 * 2.0;
//!     println!("{}", vec4); // Prints: "[8.0, 12.0, 16.0]"
//!
//!     // Common vector operations
//!     let dot_product = vec3.dot(vec4);
//!     println!("{}", dot_product); // Prints: "232"
//!
//!    let vec5 = Vec3D::new(1.0, 0.0, 0.0);
//!    let vec6 = Vec3D::new(0.0, 1.0, 0.0);
//!    let cross_product = vec5.cross(vec6);
//!    println!("{}", cross_product); // Prints: "[0.0, 0.0, 1.0]"
//!
//!    // Plus initialisations from spherical/cylindrical coordinates, rotations and more
//! }
//! ```

#[macro_use]
extern crate nearly_eq;

use std::cmp;
use std::error;
use std::f64::consts::PI;
use std::fmt;
use std::ops;

use float_cmp::ApproxEq;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Vec3DError {
    RadiusLessThanZero,
    ThetaNotWithinRange,
}

impl fmt::Display for Vec3DError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Vec3DError::RadiusLessThanZero => write!(f, "Vec3D error: given value for the radius is less than 0.0."),
            Vec3DError::ThetaNotWithinRange => write!(f, "Vec3D error: given value for theta not within the range [0, PI]"),
        }
    }
}

impl error::Error for Vec3DError {}

#[derive(Copy, Clone, Debug)]
pub struct Vec3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3D {
    pub fn new(x: f64, y: f64, z: f64) -> Vec3D {
        Vec3D { x, y, z }
    }

    /// Creates a new 3D vector from spherical coordinates
    ///
    ///  - **r** - radial distance
    ///  - **theta** - polar angle
    ///  - **phi** - azimuthal angle
    ///
    /// # Errors
    ///
    /// Will return a `Vec3DError` if:
    ///  - `r < 0`
    ///  - `theta < 0`
    ///  - `theta > PI`
    pub fn new_from_spherical_coordinates(r: f64, theta: f64, phi: f64) -> Result<Vec3D, Vec3DError> {
        if r < 0.0 {
            return Err(Vec3DError::RadiusLessThanZero);
        };

        if theta < 0.0 || theta > PI {
            return Err(Vec3DError::ThetaNotWithinRange);
        };

        let (sin_theta, cos_theta) = theta.sin_cos();
        let (sin_phi, cos_phi) = phi.sin_cos();
        let rho = r * sin_theta;

        Ok(Vec3D {
            x: rho * cos_phi,
            y: rho * sin_phi,
            z: r * cos_theta,
        })
    }

    /// Creates a new 3D vector from cylindrical coordinates
    ///
    ///  - **r** - radial distance
    ///  - **phi** - angle
    ///  - **z** - height along z-axis
    ///
    /// # Errors
    ///
    /// Will return a `Vec3DError` if:
    ///  - `r < 0`
    pub fn new_from_cylindrical_coordinates(r: f64, phi: f64, z: f64) -> Result<Vec3D, Vec3DError> {
        if r < 0.0 {
            return Err(Vec3DError::RadiusLessThanZero);
        };

        let (sin_phi, cos_phi) = phi.sin_cos();

        Ok(Vec3D {
            x: r * cos_phi,
            y: r * sin_phi,
            z,
        })
    }

    /// Returns a new vector with a value of 0.0 in each axis
    #[inline]
    pub fn zeros() -> Vec3D {
        Vec3D::new(0.0, 0.0, 0.0)
    }

    /// Returns a new vector with a value of 1.0 in each axis
    #[inline]
    pub fn ones() -> Vec3D {
        Vec3D::new(1.0, 1.0, 1.0)
    }

    /// Returns the projection of the vector in the x-axis
    #[inline]
    pub fn x_proj(&self) -> Vec3D {
        Vec3D::new(self.x, 0.0, 0.0)
    }

    /// Returns the projection of the vector in the y-axis
    #[inline]
    pub fn y_proj(&self) -> Vec3D {
        Vec3D::new(0.0, self.y, 0.0)
    }

    /// Returns the projection of the vector in the z-axis
    #[inline]
    pub fn z_proj(&self) -> Vec3D {
        Vec3D::new(0.0, 0.0, self.z)
    }

    /// Returns the vector's theta value in spherical coordinates
    #[inline]
    pub fn theta(&self) -> f64 {
        if self.x == 0.0 && self.y == 0.0 && self.z == 0.0 {
            0.0
        } else {
            (self.x * self.x + self.y * self.y).sqrt().atan2(self.z)
        }
    }

    /// Returns the vector's phi value in spherical/cylindrical coordinates
    #[inline]
    pub fn phi(&self) -> f64 {
        if self.x == 0.0 && self.y == 0.0 {
            0.0
        } else {
            self.y.atan2(self.x)
        }
    }

    /// Returns the vector's magnitude
    #[inline]
    pub fn mag(&self) -> f64 {
        self.mag2().sqrt()
    }

    /// Returns the vector's magnitude^2
    #[inline]
    pub fn mag2(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Returns a new vector of the current vector normalised
    #[inline]
    pub fn norm(&self) -> Vec3D {
        let mag2 = self.mag2();

        if mag2 == 0.0 {
            Vec3D::new(0.0, 0.0, 0.0)
        } else {
            *self * 1.0 / mag2.sqrt()
        }
    }

    /// Alias for `Vec3D::norm`
    #[inline]
    pub fn unit(&self) -> Vec3D {
        self.norm()
    }

    /// Returns the inner product of this vector with another vector
    #[inline]
    pub fn inner_product(&self, other: Vec3D) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Alias for `Vec3D::inner_product`
    #[inline]
    pub fn dot(&self, other: Vec3D) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Returns the distance between this vector and another vector
    #[inline]
    pub fn distance_to(&self, other: Vec3D) -> f64 {
        (*self - other).mag()
    }

    /// Returns the distance^2 between this vector and another vector
    #[inline]
    pub fn distance_to2(&self, other: Vec3D) -> f64 {
        (*self - other).mag2()
    }

    /// Returns the cross product of this vector with another vector
    #[inline]
    pub fn cross(&self, other: Vec3D) -> Vec3D {
        Vec3D::new(
            self.y * other.z - other.y * self.z,
            self.z * other.x - other.z * self.x,
            self.x * other.y - other.x * self.y,
        )
    }

    /// Returns the pseudo-rapidity of the vector w.r.t the z-axis
    ///
    /// See: [https://en.wikipedia.org/wiki/Pseudorapidity](https://en.wikipedia.org/wiki/Pseudorapidity)
    pub fn pseudo_rapidity(&self) -> f64 {
        let m = self.mag();

        if m.approx_eq(&0.0, 2.0 * ::std::f64::EPSILON, 2) {
            0.0
        } else if m.approx_eq(&self.z, 2.0 * ::std::f64::EPSILON, 2) {
            std::f64::INFINITY
        } else if m.approx_eq(&(-self.z), 2.0 * ::std::f64::EPSILON, 2) {
            std::f64::NEG_INFINITY
        } else {
            0.5 * ((m + self.z) / (m - self.z)).ln()
        }
    }

    /// Rotates the vector around the x-axis
    pub fn rotate_x(&mut self, angle: f64) {
        let (sinphi, cosphi) = angle.sin_cos();

        let ty = self.y * cosphi - self.z * sinphi;
        self.z = self.z * cosphi + self.y * sinphi;
        self.y = ty;
    }

    /// Rotates the vector around the y-axis
    pub fn rotate_y(&mut self, angle: f64) {
        let (sinphi, cosphi) = angle.sin_cos();

        let tz = self.z * cosphi - self.x * sinphi;
        self.x = self.x * cosphi + self.z * sinphi;
        self.z = tz;
    }

    /// Rotates the vector around the z-axis
    pub fn rotate_z(&mut self, angle: f64) {
        let (sinphi, cosphi) = angle.sin_cos();

        let tx = self.x * cosphi - self.y * sinphi;
        self.y = self.y * cosphi + self.x * sinphi;
        self.x = tx;
    }

    /// Returns true if points are approximately equal
    pub fn approx_eq(self, other: Vec3D) -> bool {
        self.x.approx_eq(&other.x, 2.0 * ::std::f64::EPSILON, 2)
            && self.y.approx_eq(&other.y, 2.0 * ::std::f64::EPSILON, 2)
            && self.z.approx_eq(&other.z, 2.0 * ::std::f64::EPSILON, 2)
    }
}

impl fmt::Display for Vec3D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{:?}, {:?}, {:?}]", self.x, self.y, self.z)
    }
}

impl cmp::PartialEq for Vec3D {
    #[inline]
    fn eq(&self, other: &Vec3D) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}
impl Eq for Vec3D {}

impl ops::Add<Vec3D> for Vec3D {
    type Output = Vec3D;

    #[inline]
    fn add(self, other: Vec3D) -> Vec3D {
        Vec3D::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl ops::AddAssign<Vec3D> for Vec3D {
    #[inline]
    fn add_assign(&mut self, other: Vec3D) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl ops::Sub<Vec3D> for Vec3D {
    type Output = Vec3D;

    #[inline]
    fn sub(self, other: Vec3D) -> Vec3D {
        Vec3D::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl ops::SubAssign<Vec3D> for Vec3D {
    #[inline]
    fn sub_assign(&mut self, other: Vec3D) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl ops::Neg for Vec3D {
    type Output = Vec3D;

    #[inline]
    fn neg(self) -> Vec3D {
        Vec3D::new(-self.x, -self.y, -self.z)
    }
}

impl ops::Mul<f64> for Vec3D {
    type Output = Vec3D;

    #[inline]
    fn mul(self, other: f64) -> Vec3D {
        Vec3D::new(self.x * other, self.y * other, self.z * other)
    }
}

impl ops::MulAssign<f64> for Vec3D {
    #[inline]
    fn mul_assign(&mut self, other: f64) {
        self.x *= other;
        self.y *= other;
        self.z *= other;
    }
}

impl ops::Div<f64> for Vec3D {
    type Output = Vec3D;

    #[inline]
    fn div(self, other: f64) -> Vec3D {
        Vec3D::new(self.x / other, self.y / other, self.z / other)
    }
}

impl ops::DivAssign<f64> for Vec3D {
    #[inline]
    fn div_assign(&mut self, other: f64) {
        self.x /= other;
        self.y /= other;
        self.z /= other;
    }
}

#[cfg(test)]
mod tests {
    use super::{Vec3D, Vec3DError};

    const PI: f64 = std::f64::consts::PI;

    #[test]
    fn create_vector() {
        // Test creation of basic vector
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        assert_eq!(vec.x, 1.0);
        assert_eq!(vec.y, 2.0);
        assert_eq!(vec.z, 3.0);
    }

    #[test]
    fn new_from_spherical_coordinates() {
        // Test creation of vector from spherical coordinates
        let vec = Vec3D::new_from_spherical_coordinates(10.0, PI / 4.0, PI / 4.0).unwrap();
        assert_nearly_eq!(vec.x, 5.0);
        assert_nearly_eq!(vec.y, 5.0);
        assert_nearly_eq!(vec.z, 7.0710678118654752);

        // Test failure condition 1
        let vec = Vec3D::new_from_spherical_coordinates(-1.0, PI / 2.0, PI / 2.0);
        assert!(vec.is_err());
        assert_eq!(vec.err().unwrap(), Vec3DError::RadiusLessThanZero);

        // Test failure condition 2
        let vec = Vec3D::new_from_spherical_coordinates(1.0, PI * 2.0, PI / 2.0);
        assert!(vec.is_err());
        assert_eq!(vec.err().unwrap(), Vec3DError::ThetaNotWithinRange);
    }

    #[test]
    fn new_from_cylindrical_coordinates() {
        // Test creation of vector from spherical coordinates
        let vec = Vec3D::new_from_cylindrical_coordinates(10.0, PI / 4.0, 20.0).unwrap();
        assert_nearly_eq!(vec.x, 7.0710678118654752);
        assert_nearly_eq!(vec.y, 7.0710678118654752);
        assert_nearly_eq!(vec.z, 20.0);

        // Test failure condition
        let vec = Vec3D::new_from_cylindrical_coordinates(-1.0, PI / 2.0, 10.0);
        assert!(vec.is_err());
        assert_eq!(vec.err().unwrap(), Vec3DError::RadiusLessThanZero);
    }

    #[test]
    fn x_proj() {
        // Test calculation of vector's projection in x-axis
        let vec = Vec3D::new(10.0, 10.0, 10.0);
        assert_eq!(vec.x_proj(), Vec3D::new(10.0, 0.0, 0.0));
    }

    #[test]
    fn y_proj() {
        // Test calculation of vector's projection in y-axis
        let vec = Vec3D::new(10.0, 10.0, 10.0);
        assert_eq!(vec.y_proj(), Vec3D::new(0.0, 10.0, 0.0));
    }

    #[test]
    fn z_proj() {
        // Test calculation of vector's projection in z-axis
        let vec = Vec3D::new(10.0, 10.0, 10.0);
        assert_eq!(vec.z_proj(), Vec3D::new(0.0, 0.0, 10.0));
    }

    #[test]
    fn theta() {
        // Test calculation of theta values from vector
        let vec = Vec3D::new(1.0, 0.0, 0.0);
        assert_eq!(vec.theta(), PI / 2.0);

        let vec = Vec3D::new(0.0, 1.0, 1.0);
        assert_eq!(vec.theta(), PI / 4.0);

        let vec = Vec3D::new(0.0, 0.0, 1.0);
        assert_eq!(vec.theta(), 0.0);

        for i in 1..8 {
            let theta = PI / i as f64;

            let vec = Vec3D::new_from_spherical_coordinates(1.0, theta, 0.0).unwrap();
            assert_nearly_eq!(vec.theta(), theta);
        }
    }

    #[test]
    fn phi() {
        // Test calculation of theta values from vector
        let vec = Vec3D::new(1.0, 0.0, 0.0);
        assert_eq!(vec.phi(), 0.0);

        let vec = Vec3D::new(0.0, 1.0, 1.0);
        assert_eq!(vec.phi(), PI / 2.0);

        let vec = Vec3D::new(0.0, 0.0, 1.0);
        assert_eq!(vec.phi(), 0.0);

        for i in 2..8 {
            let phi = (2.0 * PI) / i as f64;
            let vec = Vec3D::new_from_spherical_coordinates(1.0, PI / 2.0, phi).unwrap();
            assert_nearly_eq!(vec.phi(), phi);
        }
    }

    #[test]
    fn mag() {
        // Test calculation of vector's magnitude
        let vec = Vec3D::new(10.0, 10.0, 10.0);
        assert_nearly_eq!(vec.mag(), 17.320508075688772);
    }

    #[test]
    fn mag2() {
        // Test calculation of vector's magnitude^2
        let vec = Vec3D::new(10.0, 10.0, 10.0);
        assert_eq!(vec.mag2(), 300.0);
    }

    #[test]
    fn unit() {
        // Test calculation of vector's unit vector
        let vec = Vec3D::new(10.0, 20.0, 30.0).unit();
        assert_nearly_eq!(vec.x, 0.2672612419124243);
        assert_nearly_eq!(vec.y, 0.5345224838248487);
        assert_nearly_eq!(vec.z, 0.8017837257372731);
    }

    #[test]
    fn inner_product() {
        // Test calculation of the inner product of two vectors
        let vec1 = Vec3D::new(10.0, 20.0, 30.0);
        let vec2 = Vec3D::new(40.0, 50.0, 60.0);
        assert_eq!(vec1.inner_product(vec2), 3200.0);
        assert_eq!(vec2.inner_product(vec1), 3200.0);
    }

    #[test]
    fn dot() {
        // Test calculation of the inner product of two vectors
        let vec1 = Vec3D::new(10.0, 20.0, 30.0);
        let vec2 = Vec3D::new(40.0, 50.0, 60.0);
        assert_eq!(vec1.dot(vec2), 3200.0);
        assert_eq!(vec2.dot(vec1), 3200.0);
    }

    #[test]
    fn distance_to() {
        // Test calculation of distance between two vectors
        let vec1 = Vec3D::new(10.0, 20.0, 30.0);
        let vec2 = Vec3D::new(40.0, 50.0, 60.0);
        
        assert_nearly_eq!(vec1.distance_to(vec1), 0.0);
        assert_nearly_eq!(vec1.distance_to(vec2), 51.96152422706632);
    }

    #[test]
    fn distance_to2() {
        // Test calculation of distance^2 between two vectors
        let vec1 = Vec3D::new(10.0, 20.0, 30.0);
        let vec2 = Vec3D::new(40.0, 50.0, 60.0);
        
        assert_nearly_eq!(vec1.distance_to2(vec1), 0.0);
        assert_nearly_eq!(vec1.distance_to2(vec2), 2700.0);
    }

    #[test]
    fn cross() {
        // Test calculation of the cross product of two vectors
        let vec1 = Vec3D::new(10.0, 20.0, 30.0);
        let vec2 = Vec3D::new(40.0, 50.0, 60.0);
        assert_eq!(vec1.cross(vec2), Vec3D::new(-300.0, 600.0, -300.0));
        assert_eq!(vec2.cross(vec1), Vec3D::new(300.0, -600.0, 300.0));
    }

    #[test]
    fn pseudo_rapidity() {
        // Test when theta is PI/2, pseudo-rapidity is 0
        let vec1 = Vec3D::new(1.0, 1.0, 0.0);
        assert_eq!(vec1.pseudo_rapidity(), 0.0);

        // Test when theta is 0, pseudo-rapidity is inf
        let vec2 = Vec3D::new(0.0, 0.0, 1.0);
        assert_eq!(vec2.pseudo_rapidity(), std::f64::INFINITY);

        // Test when theta is PI/4
        let vec3 = Vec3D::new(1.0, 0.0, 1.0);
        assert_nearly_eq!(vec3.pseudo_rapidity(), 0.881373587019543025232609);
    }

    #[test]
    fn rotate_x() {
        // Test on null vector, rotate 0
        let vec = Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_x(0.0);
        assert_eq!(vec, vec2);

        // Test on null vector, rotate 180
        let vec = Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_x(PI);
        assert_eq!(vec, vec2);

        // Test on y-dir unit vector, rotate PI
        let mut vec = Vec3D::new(0.0, 1.0, 0.0);
        vec.rotate_x(PI);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, -1.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on y-dir unit vector, rotate PI/2
        let mut vec = Vec3D::new(0.0, 1.0, 0.0);
        vec.rotate_x(PI / 2.0);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 1.0);

        // Test on z-dir unit vector, rotate PI
        let mut vec = Vec3D::new(0.0, 0.0, 1.0);
        vec.rotate_x(PI);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, -1.0);

        // Test on z-dir unit vector, rotate PI/2
        let mut vec = Vec3D::new(0.0, 0.0, 1.0);
        vec.rotate_x(PI / 2.0);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, -1.0);
        assert_nearly_eq!(vec.z, 0.0);
    }

    #[test]
    fn rotate_y() {
        // Test on null vector, rotate 0
        let vec = Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_y(0.0);
        assert_eq!(vec, vec2);

        // Test on null vector, rotate 180
        let vec = Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_y(PI);
        assert_eq!(vec, vec2);

        // Test on x-dir unit vector, rotate PI
        let mut vec = Vec3D::new(1.0, 0.0, 0.0);
        vec.rotate_y(PI);
        assert_nearly_eq!(vec.x, -1.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on x-dir unit vector, rotate PI/2
        let mut vec = Vec3D::new(1.0, 0.0, 0.0);
        vec.rotate_y(PI / 2.0);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, -1.0);

        // Test on z-dir unit vector, rotate PI
        let mut vec = Vec3D::new(0.0, 0.0, 1.0);
        vec.rotate_y(PI);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, -1.0);

        // Test on z-dir unit vector, rotate PI/2
        let mut vec = Vec3D::new(0.0, 0.0, 1.0);
        vec.rotate_y(PI / 2.0);
        assert_nearly_eq!(vec.x, 1.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 0.0);
    }

    #[test]
    fn rotate_z() {
        // Test on null vector, rotate 0
        let vec = Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_z(0.0);
        assert_eq!(vec, vec2);

        // Test on null vector, rotate 180
        let vec = Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_z(PI);
        assert_eq!(vec, vec2);

        // Test on x-dir unit vector, rotate PI
        let mut vec = Vec3D::new(1.0, 0.0, 0.0);
        vec.rotate_z(PI);
        assert_nearly_eq!(vec.x, -1.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on x-dir unit vector, rotate PI/2
        let mut vec = Vec3D::new(1.0, 0.0, 0.0);
        vec.rotate_z(PI / 2.0);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 1.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on y-dir unit vector, rotate PI
        let mut vec = Vec3D::new(0.0, 1.0, 0.0);
        vec.rotate_z(PI);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, -1.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on y-dir unit vector, rotate PI/2
        let mut vec = Vec3D::new(0.0, 1.0, 0.0);
        vec.rotate_z(PI / 2.0);
        assert_nearly_eq!(vec.x, -1.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 0.0);
    }

    #[test]
    fn partial_eq() {
        let vec1 = Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = Vec3D::new(10.0, 10.0, 10.0);
        let vec3 = Vec3D::new(20.0, 20.0, 20.0);
        assert_eq!(vec1, vec2);
        assert_ne!(vec2, vec3);
    }

    #[test]
    fn add() {
        let vec1 = Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = Vec3D::new(20.0, 20.0, 20.0);
        let vec3 = Vec3D::new(30.0, 30.0, 30.0);
        assert_eq!(vec1 + vec2, vec3);
    }

    #[test]
    fn add_assign() {
        let mut vec1 = Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = Vec3D::new(20.0, 20.0, 20.0);
        let vec3 = Vec3D::new(30.0, 30.0, 30.0);

        vec1 += vec2;
        assert_eq!(vec1, vec3);
    }

    #[test]
    fn sub() {
        let vec1 = Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = Vec3D::new(20.0, 20.0, 20.0);
        let vec3 = Vec3D::new(30.0, 30.0, 30.0);
        assert_eq!(vec3 - vec2, vec1);
    }

    #[test]
    fn sub_assign() {
        let vec1 = Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = Vec3D::new(20.0, 20.0, 20.0);
        let mut vec3 = Vec3D::new(30.0, 30.0, 30.0);

        vec3 -= vec2;
        assert_eq!(vec1, vec3);
    }

    #[test]
    fn neg() {
        let vec1 = Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = Vec3D::new(-10.0, -10.0, -10.0);
        assert_eq!(-vec1, vec2);
        assert_eq!(vec1, -vec2);
    }

    #[test]
    fn mul() {
        let vec1 = Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = Vec3D::new(30.0, 30.0, 30.0);
        assert_eq!(vec1 * 3.0, vec2);
    }

    #[test]
    fn mul_assign() {
        let mut vec1 = Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = Vec3D::new(30.0, 30.0, 30.0);

        vec1 *= 3.0;
        assert_eq!(vec1, vec2);
    }

    #[test]
    fn div() {
        let vec1 = Vec3D::new(30.0, 30.0, 30.0);
        let vec2 = Vec3D::new(10.0, 10.0, 10.0);
        assert_eq!(vec1 / 3.0, vec2);
    }

    #[test]
    fn div_assign() {
        let mut vec1 = Vec3D::new(30.0, 30.0, 30.0);
        let vec2 = Vec3D::new(10.0, 10.0, 10.0);

        vec1 /= 3.0;
        assert_eq!(vec1, vec2);
    }
}
