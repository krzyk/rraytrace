use std::f64;
use std::ops::{Add, Sub, Mul, Div, Neg};

#[derive(Copy, Clone, Debug, PartialEq)]
struct Color {
    r: f64,
    g: f64,
    b: f64
}

impl Color {
    fn new(r: f64, g: f64, b: f64) -> Color {
        Color{r: r, g: g, b: b}
    }
}

impl Mul<f64> for Color {
    type Output = Color;

    fn mul(self, other: f64) -> Color {
        Color::new(self.r * other, self.g * other, self.b * other)
    }
}

impl Mul<Color> for f64 {
    type Output = Color;

    fn mul(self, other: Color) -> Color {
        Color::new(self * other.r, self * other.g, self * other.b)
    }
}

impl Add for Color {
    type Output = Color;

    fn add(self, other: Color) -> Color {
        Color::new(self.r + other.r, self.g + other.g, self.b + other.b)
    }
}

impl Sub for Color {
    type Output = Color;

    fn sub(self, other: Color) -> Color {
        Color::new(self.r - other.r, self.g - other.g, self.b - other.b)
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
struct Coords {
    x: f64,
    y: f64,
    z: f64
}

impl Coords {
    fn new(x: f64, y: f64, z: f64) -> Coords {
        Coords {x: x, y: y, z:z}
    }

    fn make_unit_vector(&self) -> Coords {
        let k = 1.0 / self.length();
        Coords::new(self.x * k, self.y * k, self.z * k)
    }

    fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    fn squared_length(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    fn dot(&self, other: Coords) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(&self, other: Coords) -> Coords {
        Coords::new(
            self.y * other.z - self.z * other.y,
            -(self.x * other.z - self.z * other.x),
            self.x * other.y - self.y * other.x
        )
    }

    fn unit_vector(self, other: Coords) -> Coords {
        self / other
    }

}

impl Add for Coords {
    type Output = Coords;

    fn add(self, other: Coords) -> Coords {
        Coords::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl Sub for Coords {
    type Output = Coords;

    fn sub(self, other: Coords) -> Coords {
        Coords::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Mul for Coords {
    type Output = Coords;

    fn mul(self, other: Coords) -> Coords {
        Coords::new(self.x * other.x, self.y * other.y, self.z * other.z)
    }
}

impl Mul<f64> for Coords {
    type Output = Coords;

    fn mul(self, other: f64) -> Coords {
        Coords::new(self.x * other, self.y * other, self.z * other)
    }
}

impl Mul<Coords> for f64 {
    type Output = Coords;

    fn mul(self, other: Coords) -> Coords {
        Coords::new(self * other.x, self * other.y, self * other.z)
    }
}

impl Div for Coords {
    type Output = Coords;

    fn div(self, other: Coords) -> Coords {
        Coords::new(self.x / other.x, self.y / other.y, self.z / other.z)
    }
}

impl Div<f64> for Coords {
    type Output = Coords;

    fn div(self, other: f64) -> Coords {
        Coords::new(self.x / other, self.y / other, self.z / other)
    }
}

impl Neg for Coords {
    type Output = Coords;

    fn neg(self) -> Coords {
        Coords::new(-self.x, -self.y, -self.z)
    }
}


#[derive(Copy, Clone, Debug, PartialEq)]
struct Ray {
    origin: Coords,
    direction: Coords
}

impl Ray {
    fn new(origin: Coords, direction: Coords) -> Ray {
        Ray {origin: origin, direction: direction}
    }

    fn point_at_parameter(self, t: f64) -> Coords {
        self.origin + t * self.direction
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
struct HitRecord {
    t: f64,
    p: Coords,
    normal: Coords
}

impl HitRecord {
    fn new() -> HitRecord {
        HitRecord {t: 0.0, p: Coords::new(0.0, 0.0, 0.0), normal: Coords::new(0.0, 0.0, 0.0)}
    }
}

trait Hittable {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64, rec: HitRecord) -> bool;
}

#[derive(Copy, Clone, Debug, PartialEq)]
struct Sphere {
    center: Coords,
    radius: f64
}

impl Sphere {

    fn new(center: Coords, radius: f64) -> Sphere {
        Sphere { center: center, radius: radius }
    }
}

impl Hittable for Sphere {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64, mut rec: HitRecord) -> bool {
        let oc = r.origin - self.center;
        let a = r.direction.dot(r.direction);
        let b = oc.dot(r.direction);
        let c = oc.dot(oc) - self.radius * self.radius;
        let discriminant = b*b - a*c;
        if discriminant > 0.0 {
            let temp = (-b - (b*b - a*c).sqrt())/a;
            if temp < t_max && temp > t_min {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = (rec.p - self.center) / self.radius;
                return true;
            }
            let other_temp = (-b + (b*b - a*c).sqrt())/a;
            if other_temp < t_max && temp > t_min {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = (rec.p - self.center) / self.radius;
                return true;
            }
        }
        return false;
    }
}

struct HitableList {
    list: Vec<Box<Hittable>>
}

impl HitableList {
    fn new(list: Vec<Box<Hittable>>) -> HitableList {
        HitableList { list: list }
    }
}

impl Hittable for HitableList {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64, mut rec: HitRecord) -> bool {
        let mut temp_rec = HitRecord::new();
        let mut hit_anything = false;
        let mut closest_so_far = t_max;
        for o in &self.list {
            if o.hit(r, t_min, closest_so_far, temp_rec) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec.clone();
            }
        }
        return hit_anything;
    }
}

fn color(r: &Ray, world: &HitableList) -> Color {
    let mut rec = HitRecord::new();
    if world.hit(r, 0.0, f64::MAX, rec) {
        0.5 * Color::new(rec.normal.x + 1.0, rec.normal.y + 1.0, rec.normal.z + 1.0)
    } else {
        let unit_direction = r.direction.make_unit_vector();
        let t = 0.5 * (unit_direction.y + 1.0);
        (1.0 - t) * Color::new(1.0, 1.0, 1.0) + t * Color::new(0.5, 0.7, 1.0)
    }
}


fn main() {
    let nx = 200;
    let ny = 100;
    println!("P3");
    println!("{} {} 255", nx, ny);
    let lower_left_corner = Coords::new(-2.0, -1.0, -1.0);
    let horizontal = Coords::new(4.0, 0.0, 0.0);
    let vertical = Coords::new(0.0, 2.0, 0.0);
    let origin = Coords::new(0.0, 0.0, 0.0);
    let world = HitableList::new(
        vec![
            Box::new(Sphere::new(Coords::new(0.0, 0.0, -1.0), 0.5)),
            Box::new(Sphere::new(Coords::new(0.0, -100.5, -1.0), 100.0)),
        ]
    );
    for j in (0..ny).rev() {
        for i in 0..nx {
            let u = i as f64 / nx as f64;
            let v = j as f64 / ny as f64;

            let r = Ray::new(
                origin,
                lower_left_corner + u * horizontal + v * vertical
            );

            let col = color(&r, &world);
            let ir = (255.99*col.r) as i32;
            let ig = (255.99*col.g) as i32;
            let ib = (255.99*col.b) as i32;
            println!("{} {} {}", ir, ig, ib);
        }
    }

    let c = Coords::new(1.0, 2.0, 3.0);
    println!("{:?}", -c);
}
