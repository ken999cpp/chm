use std::process::ExitCode;

const EPS: f64 = 1e-3;
const TOL: f64 = EPS / 2.0;
const MAX_IT: usize = 1000;

type RealFn = fn(f64) -> f64;

/// Уравнение f(x) = 0
struct Equation {
    name: &'static str,
    f: RealFn,
    df: RealFn,
    d2f: RealFn,
}

/// Отрезок [a; b]
#[derive(Clone, Copy)]
struct Seg {
    a: f64,
    b: f64,
}

fn trim_frac_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

fn non_finite(x: f64) -> &'static str {
    if x.is_nan() {
        "nan"
    } else if x > 0.0 {
        "inf"
    } else {
        "-inf"
    }
}

fn fmt_e(x: f64, prec: usize, plus: bool) -> String {
    let sign = if plus && !x.is_sign_negative() { "+" } else { "" };
    if !x.is_finite() {
        return format!("{}{}", sign, non_finite(x));
    }
    let s = format!("{:.*e}", prec, x); // Rust даёт "1.23e-4"
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    format!("{}{}e{}{:02}", sign, mant, if exp < 0 { '-' } else { '+' }, exp.abs())
}

fn fmt_g(x: f64, plus: bool) -> String {
    const P: i32 = 6;
    let sign = if plus && !x.is_sign_negative() { "+" } else { "" };
    if !x.is_finite() {
        return format!("{}{}", sign, non_finite(x));
    }
    let s = format!("{:.*e}", (P - 1) as usize, x);
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let body = if exp < -4 || exp >= P {
        format!(
            "{}e{}{:02}",
            trim_frac_zeros(mant),
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    } else {
        trim_frac_zeros(&format!("{:.*}", (P - 1 - exp) as usize, x)).to_string()
    };
    format!("{}{}", sign, body)
}

fn horner(p: &[f64], x: f64) -> Vec<f64> {
    let mut q = Vec::with_capacity(p.len());
    q.push(p[0]);
    for &c in &p[1..] {
        let prev = *q.last().unwrap();
        q.push(prev * x + c);
    }
    q
}

fn print_poly(p: &[f64]) {
    let n = p.len() - 1;
    for (i, &c) in p.iter().enumerate() {
        print!("{}*x^{} ", fmt_g(c, true), n - i);
    }
    println!();
}

fn upper_bound(p: &[f64]) -> i32 {
    for c in 0..=1000 {
        let q = horner(p, f64::from(c));
        let ok = !q.iter().any(|&v| v < 0.0);
        print!("   c = {:<3} строка Горнера:", c);
        for v in &q {
            print!(" {:9.4}", v);
        }
        println!("  {}", if ok { "<- все >= 0" } else { "" });
        if ok {
            return c;
        }
    }
    1000
}

fn poly_bounds(p: &[f64]) -> Seg {
    print!("  Верхняя граница: P(x) = ");
    print_poly(p);
    let up = upper_bound(p);
    let q: Vec<f64> = p
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            let sign = if i % 2 == 1 { -1.0 } else { 1.0 };
            c * sign + 0.0 // + 0.0 убирает -0
        })
        .collect();
    print!("  Нижняя граница: Q(x) = (-1)^n P(-x) = ");
    print_poly(&q);
    let lo = upper_bound(&q);
    println!("  => все действительные корни лежат в [{}; {}]", -lo, up);
    Seg {
        a: f64::from(-lo),
        b: f64::from(up),
    }
}

fn separate(e: &Equation, from: f64, to: f64, h: f64) -> Vec<Seg> {
    let n = ((to - from) / h).round() as i64;
    let mut r = Vec::new();
    for i in 0..n {
        let a = from + i as f64 * h;
        let b = from + (i + 1) as f64 * h;
        if (e.f)(a) * (e.f)(b) < 0.0 {
            r.push(Seg { a, b });
        }
    }
    println!(
        "  Табулирование на [{}; {}] с шагом {}: найдено отрезков с корнем: {}",
        fmt_g(from, false),
        fmt_g(to, false),
        fmt_g(h, false),
        r.len()
    );
    for s in &r {
        println!(
            "    [{:8.4}; {:8.4}]  f(a) = {:+.5}  f(b) = {:+.5}",
            s.a,
            s.b,
            (e.f)(s.a),
            (e.f)(s.b)
        );
    }
    r
}

fn narrow(e: &Equation, mut s: Seg) -> Seg {
    for _ in 0..60 {
        let need = (e.df)(s.a) * (e.df)(s.b) <= 0.0 || (e.d2f)(s.a) * (e.d2f)(s.b) <= 0.0;
        if !need {
            break;
        }
        let c = (s.a + s.b) / 2.0;
        if (e.f)(s.a) * (e.f)(c) <= 0.0 {
            s.b = c;
        } else {
            s.a = c;
        }
    }
    s
}

fn bounds1(e: &Equation, s: Seg) -> (f64, f64) {
    let x = (e.df)(s.a).abs();
    let y = (e.df)(s.b).abs();
    (x.min(y), x.max(y))
}

fn r3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0 + 0.0
}

fn check(e: &Equation, x: f64) {
    let l = (e.f)(x - EPS);
    let r = (e.f)(x + EPS);
    println!(
        "  Проверка: f({:.6}) = {};  f(x-eps) = {}, f(x+eps) = {}  => {}",
        x,
        fmt_e((e.f)(x), 2, true),
        fmt_e(l, 2, true),
        fmt_e(r, 2, true),
        if l * r <= 0.0 {
            "корень на отрезке [x-eps; x+eps], верно"
        } else {
            "ошибка"
        }
    );
}

fn iteration(e: &Equation, s: Seg) -> f64 {
    let s = narrow(e, s);
    let (m1, big_m1) = bounds1(e, s);
    let sign = if (e.df)(s.a) > 0.0 { 1.0 } else { -1.0 };
    let k = sign * big_m1;
    let q = 1.0 - m1 / big_m1; // коэф. сжатия, 0 < q < 1 для сходимости метода итерации
    println!(
        "  Метод итерации на [{:.5}; {:.5}]: phi(x) = x - f(x)/({:.5}), q = max|phi'| = {:.5} < 1",
        s.a, s.b, k, q
    );
    let mut x = s.a;
    let stop = if q > 0.0 {
        TOL * (1.0 - q) / q
    } else {
        f64::INFINITY
    };
    println!("    n  x_n          f(x_n)\n    0  {:.6}  {:+.6}", x, (e.f)(x));
    for n in 1..=MAX_IT {
        let nx = x - (e.f)(x) / k;
        println!("  {:3}  {:.6}  {:+.6}", n, nx, (e.f)(nx));
        // разность между приближениями меньше порога остановки (TOL*(1-q)/q)
        let done = (nx - x).abs() <= stop;
        x = nx;
        if done {
            break;
        }
    }
    println!(
        "  Остановка при |x_n - x_(n-1)| <= (eps/2)(1-q)/q = {}\n  Корень: x = {:.6} ~ {:.3}",
        fmt_e(stop, 2, false),
        x,
        r3(x)
    );
    check(e, x);
    x
}

fn chord(e: &Equation, s: Seg) -> f64 {
    let s = narrow(e, s);
    let (m1, _) = bounds1(e, s);
    let (c, mut x) = if (e.f)(s.a) * (e.d2f)(s.a) > 0.0 {
        (s.a, s.b)
    } else {
        (s.b, s.a)
    };
    println!(
        "  Метод хорд на [{:.5}; {:.5}]: неподвижный конец c = {:.5}, m1 = {:.5}",
        s.a, s.b, c, m1
    );
    println!(
        "    n  x_n          f(x_n)      |f|/m1\n    0  {:.6}  {:+.6}  {}",
        x,
        (e.f)(x),
        fmt_e((e.f)(x).abs() / m1, 2, false)
    );
    let mut n = 1;
    while n <= MAX_IT && (e.f)(x).abs() / m1 > TOL {
        // x_{n+1} = x_n - f(x_n) * (c - x_n) / (f(c) - f(x_n))
        x -= (e.f)(x) * (c - x) / ((e.f)(c) - (e.f)(x));
        println!(
            "  {:3}  {:.6}  {:+.6}  {}",
            n,
            x,
            (e.f)(x),
            fmt_e((e.f)(x).abs() / m1, 2, false)
        );
        n += 1;
    }
    println!("  Корень: x = {:.6} ~ {:.3}", x, r3(x));
    check(e, x);
    x
}

fn newton(e: &Equation, s: Seg) -> f64 {
    let s = narrow(e, s);
    let (m1, _) = bounds1(e, s);
    let mut x = if (e.f)(s.a) * (e.d2f)(s.a) > 0.0 { s.a } else { s.b };
    println!(
        "  Метод касательных на [{:.5}; {:.5}]: x0 = {:.5} (f*f'' > 0), m1 = {:.5}",
        s.a, s.b, x, m1
    );
    println!(
        "    n  x_n          f(x_n)      |f|/m1\n    0  {:.6}  {:+.6}  {}",
        x,
        (e.f)(x),
        fmt_e((e.f)(x).abs() / m1, 2, false)
    );
    let mut n = 1;
    while n <= MAX_IT && (e.f)(x).abs() / m1 > TOL {
        x -= (e.f)(x) / (e.df)(x); // x_{n+1} = x_n - f(x_n) / f'(x_n)
        println!(
            "  {:3}  {:.6}  {:+.6}  {}",
            n,
            x,
            (e.f)(x),
            fmt_e((e.f)(x).abs() / m1, 2, false)
        );
        n += 1;
    }
    println!("  Корень: x = {:.6} ~ {:.3}", x, r3(x));
    check(e, x);
    x
}

fn combined(e: &Equation, s: Seg) -> f64 {
    let s = narrow(e, s);
    let (mut t, mut h) = if (e.f)(s.a) * (e.d2f)(s.a) > 0.0 {
        (s.a, s.b)
    } else {
        (s.b, s.a)
    };
    println!("  Комбинированный метод на [{:.5}; {:.5}]", s.a, s.b);
    println!(
        "    n  касательная  хорда        |разность|\n    0  {:.6}    {:.6}    {}",
        t,
        h,
        fmt_e((t - h).abs(), 2, false)
    );
    let mut n = 1;
    while n <= MAX_IT && (t - h).abs() > 2.0 * TOL {
        let nh = h - (e.f)(h) * (t - h) / ((e.f)(t) - (e.f)(h));
        t -= (e.f)(t) / (e.df)(t);
        h = nh;
        println!(
            "  {:3}  {:.6}    {:.6}    {}",
            n,
            t,
            h,
            fmt_e((t - h).abs(), 2, false)
        );
        n += 1;
    }
    let x = (t + h) / 2.0;
    println!("  Корень: x = (t + h)/2 = {:.6} ~ {:.3}", x, r3(x));
    check(e, x);
    x
}

fn ft(x: f64) -> f64 {
    x * x - 2.0
}

fn dft(x: f64) -> f64 {
    2.0 * x
}

fn d2ft(_: f64) -> f64 {
    2.0
}

fn self_test() -> bool {
    let sqrt2 = 2.0_f64.sqrt();
    println!(
        "тест: f(x) = x^2 - 2 = 0, точный корень sqrt(2) = {:.6} =====",
        sqrt2
    );
    let e = Equation {
        name: "x^2 - 2 = 0",
        f: ft,
        df: dft,
        d2f: d2ft,
    };
    let s = Seg { a: 1.0, b: 2.0 };
    let r = [
        iteration(&e, s),
        chord(&e, s),
        newton(&e, s),
        combined(&e, s),
    ];
    let ok = r.iter().all(|&v| (v - sqrt2).abs() <= EPS);
    println!(
        "тест {}\n",
        if ok {
            "пройден: все 4 метода дали sqrt(2) с точностью 0.001"
        } else {
            "не пройден"
        }
    );
    ok
}

fn f21(x: f64) -> f64 {
    x * x * x - 0.1 * x * x + 0.4 * x + 2.0
}

fn df21(x: f64) -> f64 {
    3.0 * x * x - 0.2 * x + 0.4
}

fn d2f21(x: f64) -> f64 {
    6.0 * x - 0.2
}

fn f22(x: f64) -> f64 {
    (x - 1.0) * (x - 1.0) - 0.5 * x.exp()
}

fn df22(x: f64) -> f64 {
    2.0 * (x - 1.0) - 0.5 * x.exp()
}

fn d2f22(x: f64) -> f64 {
    2.0 - 0.5 * x.exp()
}

fn f3a(x: f64) -> f64 {
    x + x.log10() - 0.5
}

fn df3a(x: f64) -> f64 {
    1.0 + 1.0 / (x * 10.0_f64.ln())
}

fn d2f3a(x: f64) -> f64 {
    -1.0 / (x * x * 10.0_f64.ln())
}

fn f3b(x: f64) -> f64 {
    x * x * x + 3.0 * x + 1.0
}

fn df3b(x: f64) -> f64 {
    3.0 * x * x + 3.0
}

fn d2f3b(x: f64) -> f64 {
    6.0 * x
}

fn f3c(x: f64) -> f64 {
    x * x * x - 3.0 * x * x + 2.5
}

fn df3c(x: f64) -> f64 {
    3.0 * x * x - 6.0 * x
}

fn d2f3c(x: f64) -> f64 {
    6.0 * x - 6.0
}

fn main() -> ExitCode {
    if !self_test() {
        return ExitCode::FAILURE;
    }

    let e21 = Equation { name: "x^3 - 0.1x^2 + 0.4x + 2 = 0", f: f21, df: df21, d2f: d2f21 };
    let e22 = Equation { name: "(x-1)^2 = 0.5 e^x", f: f22, df: df22, d2f: d2f22 };
    let e3a = Equation { name: "x + lg(x) = 0.5", f: f3a, df: df3a, d2f: d2f3a };
    let e3b = Equation { name: "x^3 + 3x + 1 = 0", f: f3b, df: df3b, d2f: d2f3b };
    let e3c = Equation { name: "x^3 - 3x^2 + 2.5 = 0", f: f3c, df: df3c, d2f: d2f3c };

    println!("№ 2.1. {} =====", e21.name);
    let b = poly_bounds(&[1.0, -0.1, 0.4, 2.0]);
    for s in separate(&e21, b.a, b.b, 0.1) {
        iteration(&e21, s);
    }

    println!("\n№ 2.2. {} =====", e22.name);
    for s in separate(&e22, 0.0, 3.0, 0.1) {
        iteration(&e22, s);
    }

    println!("\n№ 3а. Метод хорд. {} =====", e3a.name);
    for s in separate(&e3a, 0.05, 2.05, 0.1) {
        chord(&e3a, s);
    }

    println!("\n№ 3б. Метод касательных. {} =====", e3b.name);
    let b3 = poly_bounds(&[1.0, 0.0, 3.0, 1.0]);
    for s in separate(&e3b, b3.a, b3.b, 0.1) {
        newton(&e3b, s);
    }

    println!("\n№ 3в. Комбинированный метод. {} =====", e3c.name);
    let b4 = poly_bounds(&[1.0, -3.0, 0.0, 2.5]);
    for s in separate(&e3c, b4.a, b4.b, 0.1) {
        combined(&e3c, s);
    }

    ExitCode::SUCCESS
}
