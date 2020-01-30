use vector3::Vector3;


pub struct RayHit {
    pub ray_length : f32,
    pub intersection_point : Vector3
}

impl RayHit {
    pub fn new(ray_length: f32, ray_origin : Vector3, ray_direction: Vector3) -> RayHit {
        let intersection_point = ray_origin + ray_direction * ray_length;

        RayHit{ray_length, intersection_point}
    }

    pub fn from_intersection(ray_length: f32, intersection_point: Vector3) -> RayHit {
        RayHit{ray_length, intersection_point}
    }
}