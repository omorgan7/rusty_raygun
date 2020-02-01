use std::vec::Vec;
use vector3::Vector3;

use traceable::Triangle;
use traceable::Sphere;

pub struct Scene {
    pub spheres: Vec<Sphere>,
    pub triangles: Vec<Triangle>,
    pub lights: Vec<Triangle>
}

const LIGHT_INTENSITY : f32 = 20.0;
pub const LIGHT_COLOUR : Vector3 = Vector3{x: LIGHT_INTENSITY, y: LIGHT_INTENSITY, z: LIGHT_INTENSITY};

impl Scene {
    pub fn new() -> Scene {

        let vertices : [[Vector3; 3]; 12] = [
            [Vector3{x: 10.0,  y: -10.0,  z: 10.0}, Vector3{x: -10.0, y: -10.0, z: 10.0},  Vector3{x: -10.0,  y: 10.0, z: 10.0}],
            [Vector3{x: 10.0,  y: -10.0,  z: 10.0}, Vector3{x: -10.0, y: 10.0, z: 10.0},  Vector3{x: 10.0,  y: 10.0, z: 10.0}],
            [Vector3{x: -10.0, y: 10.0, z: 10.0}, Vector3{x: -10.0,  y: -10.0,  z: 10.0}, Vector3{x: -10.0,  y: -10.0, z: -10.0}],
            [Vector3{x: 10.0,  y: -10.0,  z: 10.0}, Vector3{x: 10.0, y: 10.0, z: 10.0}, Vector3{x: 10.0,  y: -10.0, z: -10.0}],
            [Vector3{x: 10.0,  y: 10.0,  z: 10.0}, Vector3{x: 10.0, y: 10.0, z: -10.0}, Vector3{x: 10.0,  y: -10.0, z: -10.0}],
            [Vector3{x: -10.0, y: 10.0, z: -10.0}, Vector3{x: -10.0,  y: 10.0,  z: 10.0}, Vector3{x: -10.0,  y: -10.0, z: -10.0}],
            [Vector3{x: 10.0,  y: 10.0,  z: 10.0}, Vector3{x: -10.0, y: 10.0, z: 10.0}, Vector3{x: -10.0,  y: 10.0, z: -10.0}],
            [Vector3{x: 10.0, y: 10.0, z: -10.0}, Vector3{x: 10.0,  y: 10.0,  z: 10.0}, Vector3{x: -10.0,  y: 10.0, z: -10.0}],
            [Vector3{x: 10.0,  y: -10.0,  z: 10.0}, Vector3{x: 10.0, y: -10.0, z: -10.0}, Vector3{x: -10.0,  y: -10.0, z: -10.0}],
            [Vector3{x: 10.0, y: -10.0, z: 10.0}, Vector3{x: -10.0,  y: -10.0,  z: -10.0}, Vector3{x: -10.0,  y: -10.0, z: 10.0}],
            [Vector3{x: 10.0,  y: -10.0,  z: -10.0}, Vector3{x: -10.0,  y: 10.0, z: -10.0}, Vector3{x: -10.0, y: -10.0, z: -10.0}],
            [Vector3{x: -10.0, y: 10.0, z: -10.0}, Vector3{x: 10.0,  y: -10.0,  z: -10.0}, Vector3{x: 10.0,  y: 10.0, z: -10.0}]
        ];

        let light_verts : [[Vector3; 3]; 4] = [
            [Vector3{x: 2.0,  y: 9.0,  z: 2.0}, Vector3{x: -2.0, y: 9.0, z: 2.0}, Vector3{x: -2.0,  y: 9.0, z: -2.0}],
            [Vector3{x: -2.0, y: 9.05, z: 2.0}, Vector3{x: 2.0,  y: 9.05,  z: 2.0}, Vector3{x: -2.0,  y: 9.05, z: -2.0}],
            [Vector3{x: 2.0, y: 9.0, z: -2.0}, Vector3{x: 2.0,  y: 9.0,  z: 2.0}, Vector3{x: -2.0,  y: 9.0, z: -2.0}],
            [Vector3{x: 2.0,  y: 9.05,  z: 2.0}, Vector3{x: 2.0, y: 9.05, z: -2.0}, Vector3{x: -2.0,  y: 9.05, z: -2.0}]
        ];

        let colour = Vector3{x: 0.8, y: 0.8, z: 0.0};

        let to_tri = | v: &[Vector3; 3] | { Triangle::new(*v, colour) };

        let triangles = vertices.iter().map(to_tri).collect();
        let lights = light_verts.iter().map(to_tri).collect();

        let spheres = [
            Sphere {
                origin: Vector3{x: 0.0, y: -8.0, z: 0.0},
                colour: Vector3{x: 0.2, y: 0.8, z: 0.8},
                radius: 2.0
            },
            Sphere {
                origin: Vector3{x: 1.0, y: -4.0, z: 10.0},
                colour: Vector3{x: 0.2, y: 0.7, z: 0.1},
                radius: 6.0
            },
            Sphere {
                origin: Vector3{x: 3.0, y: -8.6, z: 2.0},
                colour: Vector3{x: 0.5, y: 0.4, z: 0.2},
                radius: 1.4
            },
            Sphere {
                origin: Vector3{x: -4.0, y: -6.8, z: -1.0},
                colour: Vector3{x: 0.1, y: 0.5, z: 0.15},
                radius: 3.2
            },
            Sphere {
                origin: Vector3{x: 5.0, y: -9.4, z: 0.5},
                colour: Vector3{x: 0.7, y: 0.2, z: 0.5},
                radius: 0.6
            }
        ];

        Scene {
            spheres: spheres.to_vec(),
            triangles: triangles,
            lights: lights
        }
    }
}