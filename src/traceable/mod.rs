use scene::Scene;
use scene::LIGHT_COLOUR;
use vector3::Vector3;
use rayhit::RayHit;

use std::f32::INFINITY;
use std::f32::consts::PI;

static MAX_ITERATIONS : i32 = 6;

pub struct Rng {
    state : u64,
    inc : u64
}

struct Basis {
    pub x : Vector3,
    pub y : Vector3,
    pub z : Vector3
}

impl Rng {
    pub fn new(seed: u64) -> Rng {
        let mut this = Rng {state: 0, inc: (0 << 1) | 1};
        this.state += seed;
        this.get_u32();
        this
    }

    pub fn get_u32(&mut self) -> u32 {
        let old_state = self.state;

        // advance internal state
        self.state = old_state * 6364136223846793005 + (self.inc | 1);

        // calculate output function (XSH RR), uses old state for max ILP
        let xor_shift : u32 = (((old_state >> 18) ^ old_state) >> 27) as u32;
        let rot : u32 = (old_state >> 59) as u32;
        let neg_rot = std::u32::MAX - rot;
        (xor_shift >> rot) | (xor_shift << (neg_rot & 31))
    }

    pub fn gen(&mut self) -> f32 {
        (self.get_u32() as f32) / (std::u32::MAX as f32)
    }
}

pub struct RandomEngine {
    pub rng: Rng
}

impl RandomEngine {
    pub fn hemispherical_ray(&mut self, normal: Basis) -> Vector3 {
        let u1 : f32 = self.rng.gen();
        let u2 : f32 = self.rng.gen();
    
        let phi = 2.0 * PI * u2;
        let sin_theta = (1.0 - u1 * u1).sqrt();

        let random = Vector3{x: sin_theta * phi.cos(), y: u1, z: sin_theta * phi.sin()};

        Vector3 {
            x: random.x * normal.x.x + random.y * normal.z.x + random.z * normal.y.x,
            y: random.x * normal.x.y + random.y * normal.z.y + random.z * normal.y.y,
            z: random.x * normal.x.z + random.y * normal.z.z + random.z * normal.y.z,
        }
    }

    pub fn cosine_ray(&mut self, normal: Basis) -> Vector3 {
        let u1 : f32 = self.rng.gen();
        let u2 : f32 = self.rng.gen();

        let phi = 2.0 * PI * u1;
        let one_minus_u2_sqrt = (1.0 - u2).sqrt();

        let random = Vector3{
            x: phi.cos() * one_minus_u2_sqrt,
            y: u2.sqrt(),
            z: phi.sin() * one_minus_u2_sqrt
        };

        Vector3 {
            x: random.x * normal.x.x + random.y * normal.z.x + random.z * normal.y.x,
            y: random.x * normal.x.y + random.y * normal.z.y + random.z * normal.y.y,
            z: random.x * normal.x.z + random.y * normal.z.z + random.z * normal.y.z,
        }
    }
}

pub trait Traceable {
    fn trace(&self,
        ray_origin: Vector3,
        ray_direction: Vector3) -> Option<RayHit>;

    fn colour(&self,
        rng : &mut RandomEngine,
        iteration : i32,
        ray_direction: Vector3,
        hit_position: Vector3,
        scene: &Scene,
        objects: &[&dyn Traceable]) -> Vector3;
}

#[derive(Clone)]
pub struct Triangle {
    pub vertices : [Vector3; 3],
    pub normal : Vector3,
    pub nx : Vector3, 
    pub ny : Vector3,
    pub edge_ba : Vector3,
    pub edge_cb : Vector3,
    pub edge_ac : Vector3,
    pub colour : Vector3,
    pub plane_offset : f32
}

#[derive(Clone)]
pub struct Sphere {
    pub origin : Vector3,
    pub colour : Vector3,
    pub radius : f32
}

trait Colourable {
   fn ambient_colour(&self) -> Vector3;
   fn normal_basis(&self, intersection_point : Vector3) -> Basis;
}

impl Colourable for Sphere {
    fn ambient_colour(&self) -> Vector3 { return self.colour;}
    fn normal_basis(&self, intersection_point : Vector3) -> Basis {

        let normal = (intersection_point - self.origin).normalised();
        let nz = normal;
        let nx = (if (nz.x).abs() > (nz.y).abs() { Vector3{x: nz.z, y: 0.0, z: -nz.x} } else { Vector3{x: 0.0, y: -nz.z, z: nz.y} }).normalised();
        let ny = nz.cross(nx);

        return Basis{x: nx, y: ny, z: nz};
    }
}
impl Colourable for Triangle {
    fn ambient_colour(&self) -> Vector3 { return self.colour;}

    fn normal_basis(&self, _: Vector3) -> Basis {
        return Basis{x: self.nx, y: self.ny, z: self.normal};
    }
}

impl Triangle {
    pub fn new(vertices : [Vector3; 3], colour : Vector3) -> Triangle {
        let a = vertices[0];
        let b = vertices[1];
        let c = vertices[2];

        let edge_ba = b - a;
        let edge_cb = c - b;
        let edge_ac = a - c;

        let normal = edge_ba.cross(-edge_ac).normalised();

        let plane_offset = normal.dot(&a);

        let nz = normal;
        let nx = (if (nz.x).abs() > (nz.y).abs() { Vector3{x: nz.z, y: 0.0, z: -nz.x} } else { Vector3{x: 0.0, y: -nz.z, z: nz.y} }).normalised();
        let ny = nz.cross(nx);

        Triangle {
            vertices: vertices,
            normal: normal,
            nx: nx,
            ny: ny,
            edge_ba: edge_ba,
            edge_cb: edge_cb,
            edge_ac: edge_ac,
            colour : colour,
            plane_offset: plane_offset
        }
    }

}

impl Traceable for Triangle {

    fn trace(&self, ray_origin: Vector3, ray_direction: Vector3) -> Option<RayHit> {
        let normal_direction = self.normal.dot(&ray_direction);

        if normal_direction > 0.0 {
            return None;
        }

        let normal_origin = self.normal.dot(&ray_origin);

        let ray_length = (self.plane_offset - normal_origin) / normal_direction;

        if ray_length < 0.0 {
            return None;
        }

        // check three edge conditions now:
        let intersection_point = ray_origin + ray_direction * ray_length;

        let edge_ba_bounds : f32 = (self.edge_ba.cross(intersection_point - self.vertices[0])).dot(&self.normal);
        if edge_ba_bounds < 0.0 {
            return None;
        }

        let edge_cb_bounds : f32 = (self.edge_cb.cross(intersection_point - self.vertices[1])).dot(&self.normal);
        if edge_cb_bounds < 0.0 {
            return None;
        }

        let edge_ac_bounds : f32 = (self.edge_ac.cross(intersection_point - self.vertices[2])).dot(&self.normal);
        if edge_ac_bounds < 0.0 {
            return None;
        }

        Some(RayHit::from_intersection(ray_length, intersection_point))
    }

    fn colour(&self, rng : &mut RandomEngine, iteration : i32, ray_direction: Vector3, hit_position: Vector3, scene: &Scene, objects: &[&dyn Traceable]) -> Vector3 {

        let surface_normal = self.normal;
        colour(self, surface_normal, rng, iteration, ray_direction, hit_position, scene, objects)
    }
}


impl Traceable for Sphere {

    fn trace(&self, ray_origin: Vector3, ray_direction: Vector3) -> Option<RayHit> {
        // Test intersection.
        // we are testing: || ray_origin - Sphere_origin ||^2 > (d . (ray_origin - Sphere_origin))^2 + r^2 

        let origin_difference = ray_origin - self.origin;
        let origin_difference_length_squared = origin_difference.length_squared();

        let projection = ray_direction.dot(&origin_difference);
        let projection_squared = projection * projection;

        let radius_squared = self.radius * self.radius;

        if origin_difference_length_squared > (projection_squared + radius_squared) {
            return None;
        }

        let second_part = (projection_squared + radius_squared - origin_difference_length_squared).sqrt();
        let ray_length = if projection < 0.0 { -projection - second_part } else { -projection + second_part };

        return if ray_length < 0.0 { None } else { Some(RayHit::new(ray_length, ray_origin, ray_direction)) };
    }

    fn colour(&self, rng : &mut RandomEngine, iteration : i32, ray_direction: Vector3, hit_position: Vector3, scene: &Scene, objects: &[&dyn Traceable]) -> Vector3 {

        let surface_normal = (hit_position - self.origin).normalised();
        // let shadow = rng.cosine_ray(ray_direction);
        colour(self, surface_normal, rng, iteration, ray_direction, hit_position, scene, objects)
    }
}

fn colour<T : Colourable>(this: &T, surface_normal : Vector3, rng : &mut RandomEngine, iteration : i32, ray_direction: Vector3, hit_position: Vector3, scene: &Scene, objects: &[&dyn Traceable]) -> Vector3 {

    if iteration == MAX_ITERATIONS {
        // our scene is a closed system
        // so this _should_ be impossible.
        // convert to assertion later on.
        return Vector3::zeroes();
    }

    // ray direction currently unused - will be used for specularity/BDRF eventually.
    let shadow = rng.cosine_ray(this.normal_basis(hit_position));

    let eps : f32 = 1e-2;

    // intersect with the nearest light
    let mut hit_length : f32 = INFINITY;
    let mut light_hit : Option<Vector3> = None;

    for light in scene.lights.iter() {
        if let Some(hit) = light.trace(hit_position, shadow) {
            if hit.ray_length + eps < hit_length {
                hit_length = hit.ray_length;

                // lights don't have different colour or intensity for now.
                light_hit = Some(LIGHT_COLOUR);
            }
        }
    }

    let mut tri_hit : Option<(&Triangle, Vector3)> = None;
    let mut sphere_hit : Option<(&Sphere, Vector3)> = None;

    for object in scene.triangles.iter() {
        let this_ptr = this as *const T as *const usize;
        let object_ptr = object as *const Triangle as *const usize;

        if this_ptr == object_ptr {
            continue;
        }

        if let Some(hit) = object.trace(hit_position, shadow) {
            if hit.ray_length + eps < hit_length {
                hit_length = hit.ray_length;

                // lights don't have different colour or intensity for now.
                tri_hit = Some((object, hit.intersection_point));
            }
        }
    }

    for object in scene.spheres.iter() {
        let this_ptr = this as *const T as *const usize;
        let object_ptr = object as *const Sphere as *const usize;

        if this_ptr == object_ptr {
            continue;
        }

        if let Some(hit) = object.trace(hit_position, shadow) {
            if hit.ray_length + eps < hit_length {
                hit_length = hit.ray_length;

                // lights don't have different colour or intensity for now.
                sphere_hit = Some((object, hit.intersection_point));
            }
        }
    }

    if let Some((object, new_hit_position)) = sphere_hit {
        let diffuse_response =  this.ambient_colour() * 0.5 * surface_normal.dot(&shadow);
        return object.colour(rng, iteration + 1, shadow, new_hit_position, scene, objects) * diffuse_response;
    }

    if let Some((object, new_hit_position)) = tri_hit {
        let diffuse_response =  this.ambient_colour() * 0.5 * surface_normal.dot(&shadow);
        return object.colour(rng, iteration + 1, shadow, new_hit_position, scene, objects) * diffuse_response;
    }

    if let Some(colour) = light_hit {
        return colour;
    }

    return Vector3::zeroes();
}
