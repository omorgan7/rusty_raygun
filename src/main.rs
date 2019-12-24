use std::fs::File;
use std::vec::Vec;
use std::io::Write;
use std::f32::consts::PI;
use std::f32::INFINITY;

extern crate rand;
use rand::Rng;
use rand::rngs::ThreadRng;

mod vector3;
use vector3::Vector3;

struct RandomEngine {
    rng: ThreadRng
}

static MAX_ITERATIONS : i32 = 4;

impl RandomEngine {
    fn hemispherical_ray(&mut self, centre: Vector3, normal: Vector3) -> Vector3 {

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
}

fn write_to_ppm(filename : &str, image : &[u8], width : i32, height : i32) -> std::io::Result<()> {

    let mut file = File::create(filename)?;

    file.write_fmt(format_args!("P6 {} {} 255 ", width, height))?;
    file.write_all(image)?;

    Ok(())
}

fn pixel_sizes(width : i32, height : i32, fov : f32, focal_length : f32) -> (f32, f32) {
    let aspect_ratio = (width as f32) / (height as f32);

    let fov_radians = fov / 360.0 * PI;

    let pixel_height = fov_radians.tan() * 2.0 * focal_length;
    let pixel_width = pixel_height * aspect_ratio;

    (pixel_width, pixel_height)
}

struct Camera
{
    up : Vector3,
    left : Vector3,
    origin : Vector3,
    image_plane_center : Vector3,
    image_plane_edge : Vector3
}

impl Camera {
    fn new(
        camera_origin : Vector3,
        look_at : Vector3,
        look_up : Vector3,
        pixel_width : f32,
        pixel_height : f32,
        focal_length : f32
    ) -> Camera {

        let forward = (look_at - camera_origin).normalised();
        
        let up = look_up.normalised();
        let left = up.cross(forward).normalised();

        let image_plane_center = camera_origin + forward * focal_length;
        let image_plane_edge = image_plane_center + left * pixel_width * 0.5 + up * pixel_height * 0.5;

        Camera {
            up: up,
            left: left,
            origin: camera_origin,
            image_plane_center: image_plane_center,
            image_plane_edge: image_plane_edge
        }
    }
}

struct RayHit {
    ray_length : f32,
    intersection_point : Vector3
}

impl RayHit {
    fn new(ray_length: f32, ray_origin : Vector3, ray_direction: Vector3) -> RayHit {
        let intersection_point = ray_origin + ray_direction * ray_length;

        RayHit{ray_length, intersection_point}
    }

    fn from_intersection(ray_length: f32, intersection_point: Vector3) -> RayHit {
        RayHit{ray_length, intersection_point}
    }
}

trait Traceable {
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
struct Sphere {
    origin : Vector3,
    colour : Vector3,
    radius : f32
}

struct Scene {
    spheres: Vec<Sphere>,
    triangles: Vec<Triangle>,
    lights: Vec<Triangle>
}

#[derive(Clone)]
struct Triangle {
    vertices : [Vector3; 3],
    normal : Vector3,
    edge_ba : Vector3,
    edge_cb : Vector3,
    edge_ac : Vector3,
    colour : Vector3,
    plane_offset : f32
}

impl Triangle {
    fn new(vertices : [Vector3; 3], colour : Vector3) -> Triangle {
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
        let shadow = rng.hemispherical_ray(hit_position, surface_normal);

        let colour = self.colour;

        let eps : f32 = 1e-2;

        // intersect with the nearest light
        let mut light_length : f32 = INFINITY;
        let mut light_hit : Option<(&dyn Traceable, Vector3)> = None;

        for light in scene.lights.iter() {
            if let Some(hit) = light.trace(hit_position, shadow) {
                if hit.ray_length + eps < light_length {
                    light_length = hit.ray_length;

                    // lights don't have different colour or intensity for now.
                    light_hit = Some((light, Vector3::ones() * 100.0));
                }
            }
        }

        let mut object_length : f32 = INFINITY;
        let mut object_hit : Option<(&dyn Traceable, Vector3)> = None;

        for object in objects.iter() {
            let self_ptr = self as *const Triangle as *const usize;
            let object_ptr = *object as *const dyn Traceable as *const usize;

            if self_ptr == object_ptr {
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

        if iteration == MAX_ITERATIONS {
            // our scene is a closed system
            // so this _should_ be impossible.
            // convert to assertion later on.
            return Vector3::zeroes();
        }

        if let Some((light, colour)) = light_hit {
            if light_length < object_length {
                return colour;
            }
        }

        match object_hit {
            None => {
                return Vector3::zeroes();
            }

            Some((object, new_hit_position)) => {
                let diffuse_response =  self.colour * 0.5 * surface_normal.dot(&shadow);
                self.colour(rng, iteration + 1, shadow, new_hit_position, scene, objects) * diffuse_response
            }
        }
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

        let surface_normal = (self.origin - hit_position).normalised();
        let shadow = rng.hemispherical_ray(hit_position, surface_normal);

        let colour = self.colour;

        let eps : f32 = 1e-2;

        // intersect with the nearest light
        let mut light_length : f32 = INFINITY;
        let mut light_hit : Option<(&dyn Traceable, Vector3)> = None;

        for light in scene.lights.iter() {
            if let Some(hit) = light.trace(hit_position, shadow) {
                if hit.ray_length + eps < light_length {
                    light_length = hit.ray_length;

                    // lights don't have different colour or intensity for now.
                    light_hit = Some((light, Vector3::ones() * 100.0));
                }
            }
        }

        let mut object_length : f32 = INFINITY;
        let mut object_hit : Option<(&dyn Traceable, Vector3)> = None;

        for object in objects.iter() {
            let self_ptr = self as *const Sphere as *const usize;
            let object_ptr = *object as *const dyn Traceable as *const usize;

            if self_ptr == object_ptr {
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

        if iteration == MAX_ITERATIONS {
            // our scene is a closed system
            // so this _should_ be impossible.
            // convert to assertion later on.
            return Vector3::zeroes();
        }

        if let Some((light, colour)) = light_hit {
            if light_length < object_length {
                return colour;
            }
        }

        match object_hit {
            None => {
                return Vector3::zeroes();
            }

            Some((object, new_hit_position)) => {
                let diffuse_response =  self.colour * 0.5 * surface_normal.dot(&shadow);
                self.colour(rng, iteration + 1, shadow, new_hit_position, scene, objects) * diffuse_response
            }
        }
    }
}

impl Scene {
    fn new() -> Scene {

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

        let light_verts : [[Vector3; 3]; 1] = [
            [Vector3{x: 2.0,  y: 9.0,  z: 2.0}, Vector3{x: -2.0, y: 9.0, z: 2.0}, Vector3{x: -2.0,  y: 9.0, z: -2.0}],
            // [Vector3{x: 2.0, y: 9.0, z: -2.0}, Vector3{x: 2.0,  y: 9.0,  z: 2.0}, Vector3{x: -2.0,  y: 9.0, z: -2.0}]
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

fn main() {
    let width = 512;
    let height = 288;
    let channels = 3;

    let field_of_view : f32 = 90.0;
    let focal_length = (width as f32) / (2.0 * (PI * field_of_view / 360.0).tan());

    let mut image : Vec<u8> = Vec::new();
    image.resize(width * height * channels as usize, 0);

    let (pixel_width, pixel_height) = pixel_sizes(width as i32, height as i32, field_of_view, focal_length);

    let world_origin = Vector3{x: 0.0, y: 0.0, z: 0.0};
    let camera_origin = Vector3{x: 0.0, y: -7.0, z: -19.9};
    let camera_look_at = world_origin;
    let camera_look_up = Vector3{x: 0.0, y: 1.0, z: 0.0};

    let camera = Camera::new(camera_origin, camera_look_at, camera_look_up, pixel_width, pixel_height, focal_length);

    let normalised_width = pixel_width / (width as f32);
    let normalised_height = pixel_height / (height as f32);

    let scene = Scene::new();

    let mut objects = Vec::new();
    scene.triangles.iter().for_each(| tri | objects.push(tri as &dyn Traceable));
    scene.spheres.iter().for_each(| s | objects.push(s as &dyn Traceable));

    let thread_count : usize = 8;

    let mut rng = RandomEngine{rng: rand::thread_rng()};

    for pixel_count in (0..(width * height * channels)).step_by(channels) {
        let x = ((pixel_count / channels) % width) as f32;
        let y = ((pixel_count / channels) / width) as f32;

        println!("Processing: {} {}", x, y);

        let ray_direction = (camera.image_plane_edge - camera.left * normalised_width * x + camera.image_plane_center - camera.up * normalised_height * y - camera.origin).normalised();

        let mut max_dist = INFINITY;
        let mut object_hit : Option<(&dyn Traceable, Vector3)> = None;

        for object in objects.iter() {
            if let Some(ray_hit) = object.trace(camera.origin, ray_direction) {
                if ray_hit.ray_length < max_dist {
                    max_dist = ray_hit.ray_length;
                    object_hit = Some((*object, ray_hit.intersection_point));
                }
            }
        }

        let mut colour = Vector3::zeroes();

        let count = 1000;

        for i in 0..count {
            if let Some((object, intersection_point)) = object_hit {
                colour = colour + object.colour(&mut rng, 0, ray_direction, intersection_point, &scene, &objects);
            }
        }

        colour = colour / (count as f32);

        let (r, g, b) = colour.to_rgb();

        image[pixel_count] = r;
        image[pixel_count + 1] = g;
        image[pixel_count + 2] = b;
    }

    write_to_ppm("out.ppm", image.as_slice(), width as i32, height as i32).expect("Failed to write image out.");
}
