#[macro_use]
extern crate nearly_eq;

use std::cmp;
use std::error;
use std::fmt;
use std::ops;

use float_cmp::ApproxEq;

const PI: f64 = std::f64::consts::PI;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Vec3DError {
    RadiusLessThanZero,
    RhoLessThanZero,
    ThetaNotWithinRange
}

impl fmt::Display for Vec3DError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Vec3DError::RadiusLessThanZero => write!(f, "Vec3D error: given value for the radius is less than 0.0."),
            Vec3DError::RhoLessThanZero => write!(f, "Vec3D error: given value for rho is less than 0.0."),
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
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Creates a new 3D vector from spherical coordinates
    ///
    /// * `r`
    /// * `theta`
    /// * `phi`
    ///
    /// Will return an Error if:
    /// * the given radius is less than 0
    /// * theta is not within the range [0, PI]
    ///
    pub fn new_from_spherical_coordinates(r: f64, theta: f64, phi: f64) -> Result<Self, Vec3DError> {
        if r < 0.0 {
            return Err(Vec3DError::RadiusLessThanZero);
        };

        if theta < 0.0 || theta > PI {
            return Err(Vec3DError::ThetaNotWithinRange);
        };

        let (sin_theta, cos_theta) = theta.sin_cos();
        let (sin_phi, cos_phi) = phi.sin_cos();
        let rho = r * sin_theta;

        Ok(Self {
            x: rho * cos_phi,
            y: rho * sin_phi,
            z: r * cos_theta,
        })
    }

    /// Creates a new 3D vector from cylindrical coordinates
    ///
    /// * `rho`
    /// * `phi`
    /// * `z`
    ///
    /// Will return an Error if:
    /// * the given rho is less than 0
    ///
    pub fn new_from_cylindrical_coordinates(rho: f64, phi: f64, z: f64) -> Result<Self, Vec3DError> {
        if rho < 0.0 {
            return Err(Vec3DError::RhoLessThanZero);
        };

        let (sin_phi, cos_phi) = phi.sin_cos();

        Ok(Self {
            x: rho * cos_phi,
            y: rho * sin_phi,
            z,
        })
    }

    /// Returns the vector's magnitude
    #[inline]
    pub fn mag(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Returns the vector's magnitude^2
    #[inline]
    pub fn mag2(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Returns the vector's unit vector
    #[inline]
    pub fn unit(&self) -> Self {
        let tot = self.mag2();

        if tot == 0.0 {
            Self::new(0.0, 0.0, 0.0)
        } else {
            *self * 1.0 / tot.sqrt()
        }
    }

    /// Returns the inner product of this vector with another
    ///
    /// * `other`
    ///
    #[inline]
    pub fn inner_product(&self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Alias for Vec3D::inner_product
    ///
    /// * `other`
    ///
    #[inline]
    pub fn dot(&self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Returns the cross product of this vector with another
    ///
    /// * `other`
    ///
    #[inline]
    pub fn cross(&self, other: Self) -> Self {
        Vec3D::new(
            self.y * other.z - other.y * self.z,
            self.z * other.x - other.z * self.x,
            self.x * other.y - other.x * self.y,
        )
    }

    /// Returns the pseudo-rapidity of the vector w.r.t the z-axis
    /// 
    /// https://en.wikipedia.org/wiki/Pseudorapidity
    /// 
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
    ///
    /// * `phi` - angle to rotate by (in radians)
    ///
    pub fn rotate_x(&mut self, phi: f64) {
        let (sinphi, cosphi) = phi.sin_cos();

        let ty = self.y * cosphi - self.z * sinphi;
        self.z = self.z * cosphi + self.y * sinphi;
        self.y = ty;
    }

    /// Rotates the vector around the y-axis
    ///
    /// * `phi` - angle to rotate by (in radians)
    ///
    pub fn rotate_y(&mut self, phi: f64) {
        let (sinphi, cosphi) = phi.sin_cos();

        let tz = self.z * cosphi - self.x * sinphi;
        self.x = self.x * cosphi + self.z * sinphi;
        self.z = tz;
    }

    /// Rotates the vector around the z-axis
    ///
    /// * `phi` - angle to rotate by (in radians)
    ///
    pub fn rotate_z(&mut self, phi: f64) {
        let (sinphi, cosphi) = phi.sin_cos();

        let tx = self.x * cosphi - self.y * sinphi;
        self.y = self.y * cosphi + self.x * sinphi;
        self.x = tx;
    }
}

impl fmt::Display for Vec3D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{:?}, {:?}, {:?}]", self.x, self.y, self.z)
    }
}

impl cmp::PartialEq for Vec3D {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}
impl Eq for Vec3D {}

impl ops::Add<Self> for Vec3D {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl ops::AddAssign<Self> for Vec3D {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl ops::Sub<Self> for Vec3D {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl ops::SubAssign<Self> for Vec3D {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl ops::Neg for Vec3D {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl ops::Mul<f64> for Vec3D {
    type Output = Self;

    #[inline]
    fn mul(self, other: f64) -> Self {
        Self::new(self.x * other, self.y * other, self.z * other)
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
    type Output = Self;

    #[inline]
    fn div(self, other: f64) -> Self {
        Self::new(self.x / other, self.y / other, self.z / other)
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
    use super::Vec3DError;

    const PI: f64 = std::f64::consts::PI;

    #[test]
    fn create_vector() {
        // Test creation of basic vector
        let vec = super::Vec3D::new(1.0, 2.0, 3.0);
        assert_eq!(vec.x, 1.0);
        assert_eq!(vec.y, 2.0);
        assert_eq!(vec.z, 3.0);
    }

    #[test]
    fn new_from_spherical_coordinates() {
        // Test creation of vector from spherical coordinates
        let vec = super::Vec3D::new_from_spherical_coordinates(10.0, PI / 4.0, PI / 4.0).unwrap();
        assert_nearly_eq!(vec.x, 5.0);
        assert_nearly_eq!(vec.y, 5.0);
        assert_nearly_eq!(vec.z, 7.0710678118654752);

        // Test failure condition 1
        let vec = super::Vec3D::new_from_spherical_coordinates(-1.0, PI/2.0, PI/2.0);
        assert!(vec.is_err());
        assert_eq!(vec.err().unwrap(), Vec3DError::RadiusLessThanZero);

        // Test failure condition 2
        let vec = super::Vec3D::new_from_spherical_coordinates(1.0, PI*2.0, PI/2.0);
        assert!(vec.is_err());
        assert_eq!(vec.err().unwrap(), Vec3DError::ThetaNotWithinRange);
    }

    #[test]
    fn new_from_cylindrical_coordinates() {
        // Test creation of vector from spherical coordinates
        let vec = super::Vec3D::new_from_cylindrical_coordinates(10.0, PI / 4.0, 20.0).unwrap();
        assert_nearly_eq!(vec.x, 7.0710678118654752);
        assert_nearly_eq!(vec.y, 7.0710678118654752);
        assert_nearly_eq!(vec.z, 20.0);

        // Test failure condition
        let vec = super::Vec3D::new_from_cylindrical_coordinates(-1.0, PI/2.0, 10.0);
        assert!(vec.is_err());
        assert_eq!(vec.err().unwrap(), Vec3DError::RhoLessThanZero);
    }

    #[test]
    fn mag() {
        // Test calculation of vector's magnitude
        let vec = super::Vec3D::new(10.0, 10.0, 10.0);
        assert_nearly_eq!(vec.mag(), 17.320508075688772);
    }

    #[test]
    fn mag2() {
        // Test calculation of vector's magnitude^2
        let vec = super::Vec3D::new(10.0, 10.0, 10.0);
        assert_eq!(vec.mag2(), 300.0);
    }

    #[test]
    fn unit() {
        // Test calculation of vector's unit vector
        let vec = super::Vec3D::new(10.0, 20.0, 30.0).unit();
        assert_nearly_eq!(vec.x, 0.2672612419124243);
        assert_nearly_eq!(vec.y, 0.5345224838248487);
        assert_nearly_eq!(vec.z, 0.8017837257372731);
    }

    #[test]
    fn inner_product() {
        // Test calculation of the inner product of two vectors
        let vec1 = super::Vec3D::new(10.0, 20.0, 30.0);
        let vec2 = super::Vec3D::new(40.0, 50.0, 60.0);
        assert_eq!(vec1.inner_product(vec2), 3200.0);
        assert_eq!(vec2.inner_product(vec1), 3200.0);
    }

    #[test]
    fn dot() {
        // Test calculation of the inner product of two vectors
        let vec1 = super::Vec3D::new(10.0, 20.0, 30.0);
        let vec2 = super::Vec3D::new(40.0, 50.0, 60.0);
        assert_eq!(vec1.dot(vec2), 3200.0);
        assert_eq!(vec2.dot(vec1), 3200.0);
    }

    #[test]
    fn cross() {
        // Test calculation of the cross product of two vectors
        let vec1 = super::Vec3D::new(10.0, 20.0, 30.0);
        let vec2 = super::Vec3D::new(40.0, 50.0, 60.0);
        assert_eq!(vec1.cross(vec2), super::Vec3D::new(-300.0, 600.0, -300.0));
        assert_eq!(vec2.cross(vec1), super::Vec3D::new(300.0, -600.0, 300.0));
    }

    #[test]
    fn pseudo_rapidity() {
        // Test when theta is PI/2, pseudo-rapidity is 0
        let vec1 = super::Vec3D::new(1.0, 1.0, 0.0);
        assert_eq!(vec1.pseudo_rapidity(), 0.0);

        // Test when theta is 0, pseudo-rapidity is inf
        let vec2 = super::Vec3D::new(0.0, 0.0, 1.0);
        assert_eq!(vec2.pseudo_rapidity(), std::f64::INFINITY);

        // Test when theta is PI/4
        let vec3 = super::Vec3D::new(1.0, 0.0, 1.0);
        assert_nearly_eq!(vec3.pseudo_rapidity(), 0.881373587019543025232609);
    }

    #[test]
    fn rotate_x() {
        // Test on null vector, rotate 0
        let vec = super::Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_x(0.0);
        assert_eq!(vec, vec2);

        // Test on null vector, rotate 180
        let vec = super::Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_x(PI);
        assert_eq!(vec, vec2);

        // Test on y-dir unit vector, rotate PI
        let mut vec = super::Vec3D::new(0.0, 1.0, 0.0);
        vec.rotate_x(PI);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, -1.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on y-dir unit vector, rotate PI/2
        let mut vec = super::Vec3D::new(0.0, 1.0, 0.0);
        vec.rotate_x(PI / 2.0);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 1.0);

        // Test on z-dir unit vector, rotate PI
        let mut vec = super::Vec3D::new(0.0, 0.0, 1.0);
        vec.rotate_x(PI);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, -1.0);

        // Test on z-dir unit vector, rotate PI/2
        let mut vec = super::Vec3D::new(0.0, 0.0, 1.0);
        vec.rotate_x(PI / 2.0);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, -1.0);
        assert_nearly_eq!(vec.z, 0.0);
    }

    #[test]
    fn rotate_y() {
        // Test on null vector, rotate 0
        let vec = super::Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_y(0.0);
        assert_eq!(vec, vec2);

        // Test on null vector, rotate 180
        let vec = super::Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_y(PI);
        assert_eq!(vec, vec2);

        // Test on x-dir unit vector, rotate PI
        let mut vec = super::Vec3D::new(1.0, 0.0, 0.0);
        vec.rotate_y(PI);
        assert_nearly_eq!(vec.x, -1.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on x-dir unit vector, rotate PI/2
        let mut vec = super::Vec3D::new(1.0, 0.0, 0.0);
        vec.rotate_y(PI / 2.0);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, -1.0);

        // Test on z-dir unit vector, rotate PI
        let mut vec = super::Vec3D::new(0.0, 0.0, 1.0);
        vec.rotate_y(PI);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, -1.0);

        // Test on z-dir unit vector, rotate PI/2
        let mut vec = super::Vec3D::new(0.0, 0.0, 1.0);
        vec.rotate_y(PI / 2.0);
        assert_nearly_eq!(vec.x, 1.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 0.0);
    }

    #[test]
    fn rotate_z() {
        // Test on null vector, rotate 0
        let vec = super::Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_z(0.0);
        assert_eq!(vec, vec2);

        // Test on null vector, rotate 180
        let vec = super::Vec3D::new(0.0, 0.0, 0.0);
        let mut vec2 = vec;
        vec2.rotate_z(PI);
        assert_eq!(vec, vec2);

        // Test on x-dir unit vector, rotate PI
        let mut vec = super::Vec3D::new(1.0, 0.0, 0.0);
        vec.rotate_z(PI);
        assert_nearly_eq!(vec.x, -1.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on x-dir unit vector, rotate PI/2
        let mut vec = super::Vec3D::new(1.0, 0.0, 0.0);
        vec.rotate_z(PI / 2.0);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, 1.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on y-dir unit vector, rotate PI
        let mut vec = super::Vec3D::new(0.0, 1.0, 0.0);
        vec.rotate_z(PI);
        assert_nearly_eq!(vec.x, 0.0);
        assert_nearly_eq!(vec.y, -1.0);
        assert_nearly_eq!(vec.z, 0.0);

        // Test on y-dir unit vector, rotate PI/2
        let mut vec = super::Vec3D::new(0.0, 1.0, 0.0);
        vec.rotate_z(PI / 2.0);
        assert_nearly_eq!(vec.x, -1.0);
        assert_nearly_eq!(vec.y, 0.0);
        assert_nearly_eq!(vec.z, 0.0);
    }

    #[test]
    fn partial_eq() {
        let vec1 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec3 = super::Vec3D::new(20.0, 20.0, 20.0);
        assert_eq!(vec1, vec2);
        assert_ne!(vec2, vec3);
    }

    #[test]
    fn add() {
        let vec1 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = super::Vec3D::new(20.0, 20.0, 20.0);
        let vec3 = super::Vec3D::new(30.0, 30.0, 30.0);
        assert_eq!(vec1 + vec2, vec3);
    }

    #[test]
    fn add_assign() {
        let mut vec1 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = super::Vec3D::new(20.0, 20.0, 20.0);
        let vec3 = super::Vec3D::new(30.0, 30.0, 30.0);

        vec1 += vec2;
        assert_eq!(vec1, vec3);
    }

    #[test]
    fn sub() {
        let vec1 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = super::Vec3D::new(20.0, 20.0, 20.0);
        let vec3 = super::Vec3D::new(30.0, 30.0, 30.0);
        assert_eq!(vec3 - vec2, vec1);
    }

    #[test]
    fn sub_assign() {
        let vec1 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = super::Vec3D::new(20.0, 20.0, 20.0);
        let mut vec3 = super::Vec3D::new(30.0, 30.0, 30.0);

        vec3 -= vec2;
        assert_eq!(vec1, vec3);
    }

    #[test]
    fn neg() {
        let vec1 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = super::Vec3D::new(-10.0, -10.0, -10.0);
        assert_eq!(-vec1, vec2);
        assert_eq!(vec1, -vec2);
    }

    #[test]
    fn mul() {
        let vec1 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = super::Vec3D::new(30.0, 30.0, 30.0);
        assert_eq!(vec1*3.0, vec2);
    }

    #[test]
    fn mul_assign() {
        let mut vec1 = super::Vec3D::new(10.0, 10.0, 10.0);
        let vec2 = super::Vec3D::new(30.0, 30.0, 30.0);

        vec1 *= 3.0;
        assert_eq!(vec1, vec2);
    }

    #[test]
    fn div() {
        let vec1 = super::Vec3D::new(30.0, 30.0, 30.0);
        let vec2 = super::Vec3D::new(10.0, 10.0, 10.0);
        assert_eq!(vec1/3.0, vec2);
    }

    #[test]
    fn div_assign() {
        let mut vec1 = super::Vec3D::new(30.0, 30.0, 30.0);
        let vec2 = super::Vec3D::new(10.0, 10.0, 10.0);

        vec1 /= 3.0;
        assert_eq!(vec1, vec2);
    }
}
