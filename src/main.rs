use std::fs::File;
use std::vec::Vec;
use std::io::Write;
use std::f32::consts::PI;
use std::f32::INFINITY;

mod vector3;
use vector3::Vector3;

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
    fn trace(&self, ray_origin: Vector3, ray_direction: Vector3) -> Option<RayHit>;
    fn colour(&self, ray_direction: Vector3, hit_position: Vector3, scene: &Scene, objects: &[&dyn Traceable]) -> Vector3;
}

#[derive(Clone)]
struct Sphere {
    origin : Vector3,
    colour : Vector3,
    radius : f32,
    diffuseness : f32,
    specularity : f32
}

struct Scene {
    sun_origin : Vector3,
    spheres: Vec<Sphere>,
    triangles: Vec<Triangle>,
}

#[derive(Clone)]
struct Triangle {
    vertices : [Vector3; 3],
    normal : Vector3,
    edge_ba : Vector3,
    edge_cb : Vector3,
    edge_ac : Vector3,
    colour : Vector3,
    plane_offset : f32,
    diffuseness : f32,
    specularity : f32
}

impl Triangle {
    fn new(vertices : [Vector3; 3], colour : Vector3, diffuseness : f32, specularity : f32) -> Triangle {
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
            plane_offset: plane_offset,
            diffuseness : diffuseness,
            specularity : specularity
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

    fn colour(&self, ray_direction: Vector3, hit_position: Vector3, scene: &Scene, objects: &[&dyn Traceable]) -> Vector3 {

        let shadow =  -(hit_position - scene.sun_origin);
        let ray_direction = shadow.normalised();
        let ray_length = shadow.length();

        let ambient_colour = self.colour * (1.0 - self.diffuseness);
        for object in objects.iter() {
            let self_ptr = self as *const Triangle as *const usize;
            let object_ptr = *object as *const dyn Traceable as *const usize;

            if self_ptr == object_ptr {
                continue;
            }

            // floating point imprecision fudge
            let eps : f32 = 1e-2;

            if let Some(hit) = object.trace(hit_position, ray_direction) {
                if hit.ray_length + eps < ray_length {
                    return ambient_colour;
                }
            }
        }

        let diffuse_response = self.normal.dot(&ray_direction);
        if (diffuse_response < 0.0) {
            return ambient_colour;
        }
        Vector3::ones() * diffuse_response * self.diffuseness + ambient_colour
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

    fn colour(&self, ray_direction: Vector3, hit_position: Vector3, scene: &Scene, objects: &[&dyn Traceable]) -> Vector3 {

        let shadow =  -(hit_position - scene.sun_origin);
        let ray_direction = shadow.normalised();
        let ray_length = shadow.length();

        for object in objects.iter() {
            let self_ptr = self as *const Sphere as *const usize;
            let object_ptr = *object as *const dyn Traceable as *const usize;

            if self_ptr == object_ptr {
                continue;
            }

            let eps : f32 = 1e-2;

            if let Some(hit) = object.trace(hit_position, ray_direction) {
                if hit.ray_length + eps < ray_length {
                    return self.ambient_colour();
                }
            }
        }

        let surface_normal = (hit_position - self.origin).normalised();
        let diffuse_response = surface_normal.dot(&ray_direction);
        if (diffuse_response < 0.0) {
            return self.ambient_colour()
        }
        Vector3::ones() * diffuse_response * self.diffuseness + self.ambient_colour()
    }
}

impl Sphere {
    fn ambient_colour(&self) -> Vector3 {
        let ambient_coefficient = 1.0 - self.diffuseness - self.specularity;

        self.colour * ambient_coefficient
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
            [Vector3{x: -10.0, y: 10.0, z: -10.0}, Vector3{x: 10.0,  y: -10.0,  z: -10.0}, Vector3{x: 10.0,  y: 10.0, z: -10.0}],
        ];

        let colour = Vector3{x: 0.8, y: 0.8, z: 0.0};

        let triangles = vertices.iter().map(| v | Triangle::new(*v, colour, 0.7, 0.1)).collect();

        let spheres = [
            Sphere {
                origin: Vector3{x: 0.0, y: -8.0, z: 0.0},
                colour: Vector3{x: 0.2, y: 0.8, z: 0.8},
                radius: 2.0,
                diffuseness: 0.7,
                specularity: 0.1,
            },
            Sphere {
                origin: Vector3{x: 1.0, y: -4.0, z: 10.0},
                colour: Vector3{x: 0.2, y: 0.7, z: 0.1},
                radius: 6.0,
                diffuseness: 0.7,
                specularity: 0.1,
            },
            Sphere {
                origin: Vector3{x: 3.0, y: -8.6, z: 2.0},
                colour: Vector3{x: 0.5, y: 0.4, z: 0.2},
                radius: 1.4,
                diffuseness: 0.7,
                specularity: 0.1,
            },
            Sphere {
                origin: Vector3{x: -4.0, y: -6.8, z: -1.0},
                colour: Vector3{x: 0.1, y: 0.5, z: 0.15},
                radius: 3.2,
                diffuseness: 0.7,
                specularity: 0.1,
            },
            Sphere {
                origin: Vector3{x: 5.0, y: -9.4, z: 0.5},
                colour: Vector3{x: 0.7, y: 0.2, z: 0.5},
                radius: 0.6,
                diffuseness: 0.7,
                specularity: 0.1,
            }
        ];

        Scene {
            sun_origin: Vector3{x: 0.0, y: 8.0, z: -7.0},
            spheres: spheres.to_vec(),
            triangles: triangles,
        }
    }
}

fn main() {
    let width = 1280;
    let height = 720;
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

    for pixel_count in (0..(width * height * channels)).step_by(channels) {
        let x = ((pixel_count / channels) % width) as f32;
        let y = ((pixel_count / channels) / width) as f32;

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
        if let Some((object, intersection_point)) = object_hit {
            colour = object.colour(ray_direction, intersection_point, &scene, &objects);
        }

        let (r, g, b) = colour.to_rgb();

        image[pixel_count] = r;
        image[pixel_count + 1] = g;
        image[pixel_count + 2] = b;
    }

    write_to_ppm("out.ppm", image.as_slice(), width as i32, height as i32).expect("Failed to write image out.");
}
