use scene::Scene;
use scene::LIGHT_COLOUR;
use vector3::Vector3;
use rayhit::RayHit;

use std::f32::INFINITY;
use std::f32::consts::PI;

static MAX_ITERATIONS : i32 = 6;

use rand::Rng;
use rand::rngs::ThreadRng;

pub struct RandomEngine {
    pub rng: ThreadRng
}

impl RandomEngine {
    pub fn hemispherical_ray(&mut self, normal: Vector3) -> Vector3 {

        let nz = normal;
        let nx = (if (nz.x).abs() > (nz.y).abs() { Vector3{x: nz.z, y: 0.0, z: -nz.x} } else { Vector3{x: 0.0, y: -nz.z, z: nz.y} }).normalised();
        let ny = nz.cross(nx);
    
        let u1 : f32 = self.rng.gen();
        let u2 : f32 = self.rng.gen();
    
        let phi = 2.0 * PI * u2;
        let sin_theta = (1.0 - u1 * u1).sqrt();

        let random = Vector3{x: sin_theta * phi.cos(), y: u1, z: sin_theta * phi.sin()};

        Vector3 {
            x: random.x * nx.x + random.y * nz.x + random.z * ny.x,
            y: random.x * nx.y + random.y * nz.y + random.z * ny.y,
            z: random.x * nx.z + random.y * nz.z + random.z * ny.z,
        }
    }

    pub fn cosine_ray(&mut self, normal: Vector3) -> Vector3 {
        let nz = normal;
        let nx = (if (nz.x).abs() > (nz.y).abs() { Vector3{x: nz.z, y: 0.0, z: -nz.x} } else { Vector3{x: 0.0, y: -nz.z, z: nz.y} }).normalised();
        let ny = nz.cross(nx);

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
            x: random.x * nx.x + random.y * nz.x + random.z * ny.x,
            y: random.x * nx.y + random.y * nz.y + random.z * ny.y,
            z: random.x * nx.z + random.y * nz.z + random.z * ny.z,
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
   fn colour(&self) -> Vector3;
}

impl Colourable for Sphere {
    fn colour(&self) -> Vector3 { return self.colour;}
}
impl Colourable for Triangle {
    fn colour(&self) -> Vector3 { return self.colour;}
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

        Triangle {
            vertices: vertices,
            normal: normal,
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
    // let shadow = rng.cosine_ray(ray_direction);
    let shadow = rng.hemispherical_ray(surface_normal);

    let eps : f32 = 1e-2;

    // intersect with the nearest light
    let mut light_length : f32 = INFINITY;
    let mut light_hit : Option<Vector3> = None;

    for light in scene.lights.iter() {
        if let Some(hit) = light.trace(hit_position, shadow) {
            if hit.ray_length + eps < light_length {
                light_length = hit.ray_length;

                // lights don't have different colour or intensity for now.
                light_hit = Some(LIGHT_COLOUR);
            }
        }
    }

    let mut object_length : f32 = INFINITY;
    let mut object_hit : Option<(&dyn Traceable, Vector3)> = None;

    for object in objects.iter() {
        let this_ptr = this as *const T as *const usize;
        let object_ptr = *object as *const dyn Traceable as *const usize;

        if this_ptr == object_ptr {
            continue;
        }

        if let Some(hit) = object.trace(hit_position, shadow) {
            if hit.ray_length + eps < object_length {
                object_length = hit.ray_length;

                // lights don't have different colour or intensity for now.
                object_hit = Some((*object, hit.intersection_point));
            }
        }
    }

    if let Some(colour) = light_hit {
        if light_length < object_length {
            return colour;
        }
    }

    match object_hit {
        None => {
            return Vector3::zeroes();
        }

        Some((object, new_hit_position)) => {
            let diffuse_response =  this.colour() * 0.5 * surface_normal.dot(&shadow);
            object.colour(rng, iteration + 1, shadow, new_hit_position, scene, objects) * diffuse_response
        }
    }
}
