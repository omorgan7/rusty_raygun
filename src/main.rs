use std::fs::File;
use std::vec::Vec;
use std::io::Write;
use std::f32::INFINITY;
use std::f32::consts::PI;


mod vector3;
use vector3::Vector3;

mod rayhit;
extern crate rand;


mod scene;
use scene::Scene;

mod traceable;
use traceable::Traceable;
use traceable::RandomEngine;

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

fn main() {
    let width = 200;
    let height = 100;
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
