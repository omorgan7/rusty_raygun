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

        let forward = (camera_origin - look_at).normalised();
        
        let up = look_up.normalised();
        let left = up.cross(&forward).normalised();

        let image_plane_center = camera_origin - forward * focal_length;
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

#[derive(Clone)]
struct Sphere {
    origin : Vector3,
    colour : Vector3,
    radius : f32,
    diffuseness : f32,
    specularity : f32
}

struct Scene
{
    sun_origin : Vector3,
    objects : [Sphere; 5]
}

impl Sphere {

    fn ambient_colour(&self) -> Vector3 {
        let ambient_coefficient = 1.0 - self.diffuseness - self.specularity;

        self.colour * ambient_coefficient
    }

    fn colour(&self, ray_direction: Vector3, hit_position: Vector3, scene: &Scene) -> Vector3 {

        let ray_direction = (scene.sun_origin - hit_position).normalised();

        for sphere in scene.objects.iter() {
            let self_ptr = self as *const Sphere;
            let sphere_ptr = sphere as *const Sphere;

            if self_ptr == sphere_ptr {
                continue;
            }

            if let Some(ray_length) = sphere.trace(hit_position, ray_direction) {
                return self.ambient_colour();
            }
        }

        let surface_normal = (hit_position - self.origin).normalised();
        let diffuse_response = surface_normal.dot(&ray_direction);
        
        Vector3::ones() * diffuse_response + self.ambient_colour()
    }

    fn trace(&self, ray_origin: Vector3, ray_direction: Vector3) -> Option<f32> {
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

        return if ray_length < 0.0 { None } else { Some(ray_length) };
    }
}

impl Scene {
    fn new() -> Scene {
        let spheres = [
            Sphere {
                origin: Vector3{x: 0.0, y: -3.0, z: 0.0},
                colour: Vector3{x: 0.8, y: 0.8, z: 0.0},
                radius: 1.0,
                diffuseness: 0.7,
                specularity: 0.1,
            },
            Sphere {
                origin: Vector3{x: 1.0, y: 1.0, z: 10.0},
                colour: Vector3{x: 0.2, y: 0.7, z: 0.1},
                radius: 6.0,
                diffuseness: 0.7,
                specularity: 0.1,
            },
            Sphere {
                origin: Vector3{x: 3.0, y: -4.0, z: 2.0},
                colour: Vector3{x: 0.5, y: 0.4, z: 0.2},
                radius: 1.4,
                diffuseness: 0.7,
                specularity: 0.1,
            },
            Sphere {
                origin: Vector3{x: -4.0, y: 6.0, z: -1.0},
                colour: Vector3{x: 0.1, y: 0.5, z: 0.15},
                radius: 3.2,
                diffuseness: 0.7,
                specularity: 0.1,
            },
            Sphere {
                origin: Vector3{x: 5.0, y: 1.0, z: 0.5},
                colour: Vector3{x: 0.7, y: 0.2, z: 0.5},
                radius: 0.6,
                diffuseness: 0.7,
                specularity: 0.1,
            }
        ];

        Scene{sun_origin: Vector3{x: 2.0, y: -2.0, z: -4.0}, objects: spheres}
    }
}

fn main() {
    let width = 1280;
    let height = 720;
    let channels = 3;

    let field_of_view : f32 = 90.0;
    let focal_length : f32 = 10.0;

    let mut image : Vec<u8> = Vec::new();
    image.resize(width * height * channels as usize, 0);

    let (pixel_width, pixel_height) = pixel_sizes(width as i32, height as i32, field_of_view, focal_length);

    let world_origin = Vector3{x: 0.0, y: 0.0, z: 0.0};
    let camera_origin = Vector3{x: 0.0, y: 0.0, z: -10.0};
    let camera_look_at = world_origin;
    let camera_look_up = Vector3{x: 0.0, y: 1.0, z: 0.0};

    let camera = Camera::new(camera_origin, camera_look_at, camera_look_up, pixel_width, pixel_height, focal_length);

    let normalised_width = pixel_width / (width as f32);
    let normalised_height = pixel_height / (height as f32);

    let scene = Scene::new();


    for pixel_count in (0..(width * height * channels)).step_by(channels) {
        let x = ((pixel_count / channels) % width) as f32;
        let y = ((pixel_count / channels) / width) as f32;

        let ray_direction = (camera.image_plane_edge - camera.left * normalised_width * x + camera.image_plane_center - camera.up * normalised_height * y - camera.origin).normalised();

        let mut max_dist = INFINITY;
        let mut object_hit : Option<(&Sphere, f32)> = None;


        for sphere in scene.objects.iter() {
            if let Some(ray_length) = sphere.trace(camera.origin, ray_direction) {

                if ray_length < max_dist {
                    max_dist = ray_length;
                    object_hit = Some((&sphere, ray_length));
                }
            }
        }

        let mut colour = Vector3::zeroes();
        if let Some((sphere, ray_length)) = object_hit {
            
            let intersection_point = camera.origin + ray_direction * ray_length;
            colour = sphere.colour(ray_direction, intersection_point, &scene);
        }

        let (r, g, b) = colour.to_rgb();

        image[pixel_count] = r;
        image[pixel_count + 1] = g;
        image[pixel_count + 2] = b;
    }

    write_to_ppm("out.ppm", image.as_slice(), width as i32, height as i32).expect("Failed to write image out.");
}
