type Point = [f64; 2];
type Jacobian = [[f64; 2]; 2];

const EPS: f64 = 0.001;
const TOL: f64 = EPS / 100.0;
const MAX_STEPS: usize = 10_000;

#[derive(Clone, Copy)]
struct Solution {
    point: Point,
    steps: usize,
}

fn norm_inf([x, y]: Point) -> f64 {
    x.abs().max(y.abs())
}

fn dot([x, y]: Point, [u, v]: Point) -> f64 {
    x * u + y * v
}

// Одновременное обновление компонент и оценка Банаха в норме max.
fn fixed_point(
    mut point: Point,
    q: f64,
    map: impl Fn(Point) -> Point,
) -> Result<(Solution, f64), &'static str> {
    for steps in 1..=MAX_STEPS {
        let next = map(point);
        let difference = [next[0] - point[0], next[1] - point[1]];
        let bound = q / (1.0 - q) * norm_inf(difference);
        point = next;
        if bound <= TOL {
            return Ok((Solution { point, steps }, bound));
        }
    }
    Err("Превышен предел числа простых итераций.")
}

// Градиентный спуск для S = (f1² + f2²)/2 с условием Армихо.
fn gradient_descent(
    mut point: Point,
    evaluate: impl Fn(Point) -> (Point, Jacobian),
) -> Result<Solution, &'static str> {
    for steps in 0..=MAX_STEPS {
        let (f, jac) = evaluate(point);
        if norm_inf(f) <= TOL {
            return Ok(Solution { point, steps });
        }

        let s = 0.5 * dot(f, f);
        let grad = [
            jac[0][0] * f[0] + jac[1][0] * f[1],
            jac[0][1] * f[0] + jac[1][1] * f[1],
        ];
        let grad_sq = dot(grad, grad);
        if !grad_sq.is_finite() || grad_sq == 0.0 {
            return Err("Спуск остановился при ненулевой невязке.");
        }

        let mut alpha = 1.0;
        loop {
            let trial = [point[0] - alpha * grad[0], point[1] - alpha * grad[1]];
            let (new_f, _) = evaluate(trial);
            let new_s = 0.5 * dot(new_f, new_f);
            if new_s <= s - 1e-4 * alpha * grad_sq && new_s < s {
                point = trial;
                break;
            }
            alpha *= 0.5;
            if alpha < 2.0_f64.powi(-40) {
                return Err("Не удалось подобрать шаг градиентного спуска.");
            }
        }
    }
    Err("Превышен предел числа шагов градиентного спуска.")
}

fn first_map([x, y]: Point) -> Point {
    [-0.4 - y.sin(), 0.5 * (x + 1.0).cos()]
}

fn first_residual([x, y]: Point) -> Point {
    [2.0 * y - (x + 1.0).cos(), x + y.sin() + 0.4]
}

fn second_system([x, y]: Point) -> (Point, Jacobian) {
    let c = (x + y).cos();
    (
        [(x + y).sin() - 1.5 * x - 0.1, x * x + y * y - 1.0],
        [[c - 1.5, c], [2.0 * x, 2.0 * y]],
    )
}

fn on_circle(x: f64, sign: f64) -> f64 {
    let y = sign * (1.0 - x * x).sqrt();
    second_system([x, y]).0[0]
}

// На нужных полуинтервалах |g'(x)| >= 0.5 и |y_c'(x)| <= 1.1.
fn gradient_error([x, y]: Point, [f1, f2]: Point, sign: f64) -> Option<f64> {
    let valid_x = if sign > 0.0 {
        (0.0..=0.6).contains(&x)
    } else {
        (-11.0 / 15.0..=0.0).contains(&x)
    };
    if !valid_x || sign * y <= 0.0 {
        return None;
    }
    let circle_error = 1.5 * f2.abs();
    let x_error = 2.0 * (f1.abs() + circle_error);
    Some(x_error.max(circle_error + 1.1 * x_error))
}

fn print_result(solution: Solution, residual: Point, error: f64) -> Result<(), &'static str> {
    if !error.is_finite() || error > EPS || norm_inf(residual) > EPS {
        return Err("Точность решения не подтверждена.");
    }
    let [x, y] = solution.point;
    println!("  Шагов: {}", solution.steps);
    println!("  x = {x:.8}; y = {y:.8}");
    println!("  F = ({:+.4e}; {:+.4e})", residual[0], residual[1]);
    println!("  ||F||_inf = {:.4e}", norm_inf(residual));
    println!("  Верхняя граница ошибки координат = {error:.4e}");
    Ok(())
}

fn self_test() -> Result<(), &'static str> {
    let exact = [1.0, -2.0];
    let (by_iteration, _) =
        fixed_point([0.0, 0.0], 0.5, |[x, y]| [(x + 1.0) / 2.0, (y - 2.0) / 2.0])?;
    let by_gradient = gradient_descent([0.0, 0.0], |[x, y]| {
        ([x - 1.0, 2.0 * (y + 2.0)], [[1.0, 0.0], [0.0, 2.0]])
    })?;
    println!("Тест с точным решением (1; -2):");
    for (name, result) in [("Итерации", by_iteration), ("Градиент", by_gradient)] {
        let error = norm_inf([result.point[0] - exact[0], result.point[1] - exact[1]]);
        if error > EPS {
            return Err("Контрольный пример не пройден.");
        }
        println!("  {name}: шагов {}, ошибка {error:.4e}", result.steps);
    }
    Ok(())
}

fn main() -> Result<(), &'static str> {
    println!("Лабораторная работа 3, вариант 8. eps = {EPS}");
    self_test()?;

    println!("\nв) Метод простых итераций");
    let start = [-0.75, 0.45];
    let q = 0.4_f64.cos();
    println!("  Старт: ({:.2}; {:.2}), q = {q:.8}", start[0], start[1]);
    let (solution, bound) = fixed_point(start, q, first_map)?;
    print_result(solution, first_residual(solution.point), bound)?;

    println!("\nг) Градиентный метод");
    for (sign, a, b) in [(1.0, 0.5, 0.6), (-1.0, -11.0 / 15.0, -0.7)] {
        let left = on_circle(a, sign);
        let right = on_circle(b, sign);
        if !(left > 0.0 && right < 0.0) {
            return Err("На дуге не подтверждена смена знака функции.");
        }
        let x0 = (a + b) / 2.0;
        let start = [x0, sign * (1.0_f64 - x0 * x0).sqrt()];
        println!("  Ветвь y со знаком {sign:+.0}: [{a:.6}; {b:.6}]");
        println!("  g(a) = {left:+.6}; g(b) = {right:+.6}");
        println!("  Старт: ({:.6}; {:.6})", start[0], start[1]);
        let solution = gradient_descent(start, second_system)?;
        let residual = second_system(solution.point).0;
        let bound = gradient_error(solution.point, residual, sign)
            .ok_or("Приближение вышло за область оценки точности.")?;
        print_result(solution, residual, bound)?;
    }
    Ok(())
}
