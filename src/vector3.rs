extern crate num;

use std::ops;
use self::num::clamp;

#[derive(Clone, Copy, Debug)]
pub struct Vector3
{
    pub x : f32,
    pub y : f32,
    pub z : f32
}

impl Vector3 {
    pub fn zeroes() -> Vector3 {
        Vector3{x: 0.0, y: 0.0, z: 0.0}
    }

    pub fn ones() -> Vector3 {
        Vector3{x: 1.0, y: 1.0, z: 1.0}
    }

    pub fn dot(&self, with : &Vector3) -> f32 {
        self.x * with.x + self.y * with.y + self.z * with.z
    }

    pub fn cross(&self, with : &Vector3) -> Vector3 {
        let first  = self.y * with.z - self.z * with.y;
        let second = self.z * with.x - self.x * with.z;
        let third  = self.x * with.y - self.y * with.x;

        Vector3{x : first, y : second, z : third}
    }

    pub fn square_distance(&self, to : &Vector3) -> f32 {
        let x_diff = self.x - to.x;
        let y_diff = self.y - to.y;
        let z_diff = self.z - to.z;

        x_diff * x_diff + y_diff * y_diff + z_diff * z_diff
    }

    pub fn length_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn length(&self) -> f32 {
        self.length_squared().sqrt()
    }

    pub fn normalise(&mut self) {
        let vector_length = self.length();
        
        self.x /= vector_length;
        self.y /= vector_length;
        self.z /= vector_length;
    }

    pub fn normalised(&self) -> Vector3 {
        let mut tmp = self.clone();
        tmp.normalise();

        tmp
    }

    pub fn to_rgb(&self) -> (u8, u8, u8) {
        // the 1.0 == 255 encoding shall be used.
        let divisor : f32 = 255.0;
        let normalised = self * divisor;

        let x_clamped = clamp(normalised.x, 0.0, 255.0) as u8;
        let y_clamped = clamp(normalised.y, 0.0, 255.0) as u8;
        let z_clamped = clamp(normalised.z, 0.0, 255.0) as u8;

        (x_clamped, y_clamped, z_clamped)
    }
}

impl ops::Add<Vector3> for Vector3 {
    type Output = Vector3;

    fn add(self, other: Vector3) -> Vector3 {
        Vector3{x: self.x + other.x, y: self.y + other.y, z: self.z + other.z}
    }
}

impl<'a> ops::Add<Vector3> for &'a Vector3 {
    type Output = Vector3;

    fn add(self, other: Vector3) -> Vector3 {
        Vector3{x: self.x + other.x, y: self.y + other.y, z: self.z + other.z}
    }
}

impl ops::Sub<Vector3> for Vector3 {
    type Output = Vector3;

    fn sub(self, other: Vector3) -> Vector3 {
        Vector3{x: self.x - other.x, y: self.y - other.y, z: self.z - other.z}
    }
}

impl<'a> ops::Sub<Vector3> for &'a Vector3 {
    type Output = Vector3;

    fn sub(self, other: Vector3) -> Vector3 {
        Vector3{x: self.x - other.x, y: self.y - other.y, z: self.z - other.z}
    }
}

impl ops::Mul<Vector3> for Vector3 {
    type Output = Vector3;

    fn mul(self, other: Vector3) -> Vector3 {
        Vector3{x: self.x * other.x, y: self.y * other.y, z: self.z * other.z}
    }
}

impl<'a> ops::Mul<Vector3> for &'a Vector3 {
    type Output = Vector3;

    fn mul(self, other: Vector3) -> Vector3 {
        Vector3{x: self.x * other.x, y: self.y * other.y, z: self.z * other.z}
    }
}

impl ops::Div<Vector3> for Vector3 {
    type Output = Vector3;

    fn div(self, other: Vector3) -> Vector3 {
        Vector3{x: self.x / other.x, y: self.y / other.y, z: self.z / other.z}
    }
}

impl<'a> ops::Div<Vector3> for &'a Vector3 {
    type Output = Vector3;

    fn div(self, other: Vector3) -> Vector3 {
        Vector3{x: self.x / other.x, y: self.y / other.y, z: self.z / other.z}
    }
}

impl ops::Mul<f32> for Vector3 {
    type Output = Vector3;

    fn mul(self, multiplier: f32) -> Vector3 {
        Vector3{x: self.x * multiplier, y: self.y * multiplier, z: self.z * multiplier}
    }
}

impl<'a> ops::Mul<f32> for &'a Vector3 {
    type Output = Vector3;

    fn mul(self, multiplier: f32) -> Vector3 {
        Vector3{x: self.x * multiplier, y: self.y * multiplier, z: self.z * multiplier}
    }
}

impl ops::Div<f32> for Vector3 {
    type Output = Vector3;

    fn div(self, divisor: f32) -> Vector3 {
        Vector3{x: self.x / divisor, y: self.y / divisor, z: self.z / divisor}
    }
}

impl<'a> ops::Div<f32> for &'a Vector3 {
    type Output = Vector3;

    fn div(self, divisor: f32) -> Vector3 {
        Vector3{x: self.x / divisor, y: self.y / divisor, z: self.z / divisor}
    }
}

impl ops::Neg for Vector3 {
    type Output = Vector3;

    fn neg(self) -> Vector3 {
        Vector3{x: -self.x, y: -self.y, z: -self.z}
    }
}