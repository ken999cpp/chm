use std::hint::black_box;
use std::time::Instant;

const N: usize = 4;
type Matrix = [[f64; N]; N];
type Vector = [f64; N];

const A: Matrix = [
    [5.0, 1.0, -1.0, 2.0],
    [0.5, -4.0, 2.0, -1.0],
    [-1.0, 0.2, -6.0, 3.0],
    [1.0, -0.5, 3.0, -6.0],
];
const B: Vector = [30.5, -20.0, 6.5, 1.2];
const EPS: f64 = 1e-8; // Требуемая абсолютная погрешность в норме max.
const MAX_ITER: usize = 10_000;
const REPEATS: usize = 100_000;

// Норма max: максимальный модуль компоненты вектора.
fn norm_inf(v: &Vector) -> f64 {
    v.iter().map(|value| value.abs()).fold(0.0, f64::max)
}

fn distance(x: &Vector, y: &Vector) -> f64 {
    let mut difference = [0.0; N];
    for i in 0..N {
        difference[i] = x[i] - y[i];
    }
    norm_inf(&difference)
}

fn mat_vec(a: &Matrix, x: &Vector) -> Vector {
    let mut result = [0.0; N];
    for i in 0..N {
        for j in 0..N {
            result[i] += a[i][j] * x[j];
        }
    }
    result
}

// Невязка r = Ax - b: проверка решения подстановкой в исходную систему.
fn residual(a: &Matrix, b: &Vector, x: &Vector) -> Vector {
    let mut result = mat_vec(a, x);
    for i in 0..N {
        result[i] -= b[i];
    }
    result
}

// Метод Халецкого: A = LU, L нижняя треугольная, диагональ U равна 1.
// Перестановки строк для матрицы этого варианта не требуются.
fn khaletsky(a: &Matrix, b: &Vector) -> Result<(Vector, Matrix, Matrix), String> {
    let mut l = [[0.0; N]; N];
    let mut u = [[0.0; N]; N];

    for j in 0..N {
        u[j][j] = 1.0;

        // Сначала вычисляем j-й столбец L, включая диагональный элемент.
        for i in j..N {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[i][k] * u[k][j];
            }
            l[i][j] = a[i][j] - sum;
        }

        if l[j][j].abs() < 1e-14 {
            return Err(format!(
                "Малый ведущий элемент L[{}][{}]: схема без перестановок неприменима.",
                j + 1,
                j + 1
            ));
        }

        // Затем вычисляем элементы j-й строки U правее диагонали.
        for i in (j + 1)..N {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[j][k] * u[k][i];
            }
            u[j][i] = (a[j][i] - sum) / l[j][j];
        }
    }

    // Прямой ход: Ly = b, неизвестные находим сверху вниз.
    let mut y = [0.0; N];
    for i in 0..N {
        let mut sum = 0.0;
        for j in 0..i {
            sum += l[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / l[i][i];
    }

    // Обратный ход: Ux = y, неизвестные находим снизу вверх.
    let mut x = [0.0; N];
    for i in (0..N).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..N {
            sum += u[i][j] * x[j];
        }
        x[i] = y[i] - sum; // Делить на u[i][i] не нужно: это 1.
    }
    Ok((x, l, u))
}

// Строим эквивалентную систему x = alpha*x + beta.
fn iteration_form(a: &Matrix, b: &Vector) -> Result<(Matrix, Vector, f64), String> {
    let mut alpha = [[0.0; N]; N];
    let mut beta = [0.0; N];
    let mut q: f64 = 0.0;
    for i in 0..N {
        if a[i][i].abs() < 1e-14 {
            return Err("Нулевой элемент на диагонали: нужно преобразовать систему.".into());
        }
        beta[i] = b[i] / a[i][i];
        let mut row_sum = 0.0;
        for j in 0..N {
            if i != j {
                alpha[i][j] = -a[i][j] / a[i][i];
                row_sum += alpha[i][j].abs();
            }
        }
        q = q.max(row_sum);
    }
    if q >= 1.0 {
        return Err(format!(
            "q = {q}: достаточное условие сходимости q < 1 не выполнено."
        ));
    }
    Ok((alpha, beta, q))
}

// Простые итерации (Якоби): все новые компоненты зависят только от старого x.
fn simple_iterations(a: &Matrix, b: &Vector, eps: f64) -> Result<(Vector, usize, f64), String> {
    if !eps.is_finite() || eps <= 0.0 {
        return Err("Точность должна быть конечной и положительной.".into());
    }
    let (alpha, beta, q) = iteration_form(a, b)?;
    let mut x = [0.0; N]; // Начальное приближение x^(0).

    for iteration in 1..=MAX_ITER {
        let mut next = beta;
        for i in 0..N {
            for j in 0..N {
                next[i] += alpha[i][j] * x[j];
            }
        }
        if next.iter().any(|value| !value.is_finite()) {
            return Err("Переполнение при вычислении итераций.".into());
        }

        // Апостериорная оценка ||x* - next|| <= q/(1-q) * ||next-x||.
        let error_bound = q / (1.0 - q) * distance(&next, &x);
        x = next;
        if error_bound <= eps {
            return Ok((x, iteration, error_bound));
        }
    }
    Err(format!("За {MAX_ITER} итераций точность не достигнута."))
}

fn print_vector(name: &str, x: &Vector) {
    println!(
        "{name} = [{:.12}, {:.12}, {:.12}, {:.12}]",
        x[0], x[1], x[2], x[3]
    );
}

fn print_matrix(name: &str, a: &Matrix) {
    println!("{name}:");
    for row in a {
        println!(
            "  [{:15.10}, {:15.10}, {:15.10}, {:15.10}]",
            row[0], row[1], row[2], row[3]
        );
    }
}

// Среднее время одного полного решения, без вывода в консоль.
// black_box препятствует удалению повторных расчетов оптимизатором.
fn average_ns(mut solve: impl FnMut() -> Vector) -> f64 {
    for _ in 0..100 {
        black_box(solve()); // Небольшой прогрев.
    }
    let start = Instant::now();
    for _ in 0..REPEATS {
        black_box(solve());
    }
    start.elapsed().as_secs_f64() * 1e9 / REPEATS as f64
}

fn main() -> Result<(), String> {
    println!("Лабораторная работа 2, вариант 8. Точность eps = {EPS:e}");
    print_matrix("A", &A);
    print_vector("b", &B);

    println!("\nПроверка строгого диагонального преобладания:");
    let mut dominant = true;
    for i in 0..N {
        let off_diagonal: f64 = (0..N).filter(|&j| j != i).map(|j| A[i][j].abs()).sum();
        println!("Строка {}: {} > {}", i + 1, A[i][i].abs(), off_diagonal);
        dominant &= A[i][i].abs() > off_diagonal;
    }
    if !dominant {
        return Err("Строгое диагональное преобладание не подтверждено.".into());
    }
    println!("Матрица невырожденная; система имеет единственное решение.");

    let (direct, l, u) = khaletsky(&A, &B)?;
    print_matrix("\nL", &l);
    print_matrix("U", &u);
    let determinant: f64 = (0..N).map(|i| l[i][i]).product();
    println!("det(A) = {determinant:.12}");
    print_vector("Халецкий", &direct);
    let direct_residual = residual(&A, &B, &direct);
    print_vector("Ax - b", &direct_residual);
    println!("||Ax - b||_inf = {:.6e}", norm_inf(&direct_residual));

    let (alpha, beta, q) = iteration_form(&A, &B)?;
    print_matrix("\nalpha", &alpha);
    print_vector("beta", &beta);
    println!("q = ||alpha||_inf = {q:.12} < 1: итерации сходятся.");
    let (iterative, iterations, error_bound) = simple_iterations(&A, &B, EPS)?;
    print_vector("Простые итерации", &iterative);
    println!("Количество итераций: {iterations}");
    println!("Оценка абсолютной погрешности: {error_bound:.6e}");
    println!(
        "||Ax - b||_inf = {:.6e}",
        norm_inf(&residual(&A, &B, &iterative))
    );
    println!(
        "Различие методов ||x_итерации - x_Халецкий||_inf = {:.6e}",
        distance(&iterative, &direct)
    );

    // В обоих замерах учитываем подготовку метода и вычисление решения.
    let direct_ns = average_ns(|| khaletsky(black_box(&A), black_box(&B)).unwrap().0);
    let iterative_ns = average_ns(|| {
        simple_iterations(black_box(&A), black_box(&B), black_box(EPS))
            .unwrap()
            .0
    });
    println!("\nСреднее по {REPEATS} полным решениям, без печати:");
    println!("Халецкий: {direct_ns:.2} нс");
    println!("Итерации: {iterative_ns:.2} нс");
    println!(
        "Отношение времени итераций к Халецкому: {:.2}",
        iterative_ns / direct_ns
    );
    Ok(())
}