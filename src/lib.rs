//! calculate Sturm chain and find all real roots by Sturm method.
//!
//! ```
//! use find_real_roots_of_polynomial::SturmChain;
//! use polynomial_ring::{polynomial, Polynomial};
//! use num::{BigRational, BigInt};
//! let v = [-2, 0, 1].iter().map(|x| BigRational::from(BigInt::from(*x))).collect::<Vec<_>>();
//! let f = Polynomial::new(v); // f=x^2-2=(x-√2)(x+√2)
//! let limit = BigRational::new(BigInt::from(1), BigInt::from(16)); // 1/16
//! let sc = SturmChain::<BigRational>::new(f);
//! let roots = sc.find_all_real_roots(&limit);
//! assert_eq!(roots.len(), 2);
//! assert!(&roots[0].1 - &roots[0].0 < limit);
//! assert!(&roots[1].1 - &roots[1].0 < limit);
//! // -√2 ∈ [-93/64, -45/32] = [-1.453125, -1.40625]
//! //  √2 ∈ [45/32, 93/64] = [1.40625, 1.453125]
//! assert_eq!(roots[0].0, BigRational::new(BigInt::from(-93), BigInt::from(64)));
//! assert_eq!(roots[0].1, BigRational::new(BigInt::from(-45), BigInt::from(32)));
//! assert_eq!(roots[1].0, BigRational::new(BigInt::from(45), BigInt::from(32)));
//! assert_eq!(roots[1].1, BigRational::new(BigInt::from(93), BigInt::from(64)));
//! ```
use num_traits::{One, Zero};
use polynomial_ring::Polynomial;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

fn abs<K>(x: K) -> K
where
    K: Sized + Ord + Zero + Neg<Output = K>,
{
    if x > K::zero() {
        x
    } else {
        -x
    }
}

fn keep_sign_monic<K>(mut f: Polynomial<K>) -> Polynomial<K>
where
    K: Sized + Clone + Ord + Zero + Neg<Output = K> + for<'x> DivAssign<&'x K>,
{
    if let Some(c) = f.lc() {
        let c = abs(c.to_owned());
        f /= c;
    }
    f
}

fn get_range<K>(mut f: Polynomial<K>) -> K
where
    K: Sized + Clone + Ord + Zero + One + Neg<Output = K> + for<'x> DivAssign<&'x K>,
{
    f.monic();
    f.coeffs().iter().cloned().map(abs).max().unwrap() + K::one()
}

/// Find root of `f` from `(left, right]` by bisection method.
///
/// Require: `f` has exact one solution on `(left, right]`.
/// Iteration stops if `right - left < limit`.
/// ```
/// use find_real_roots_of_polynomial::bisection_method;
/// use polynomial_ring::{polynomial, Polynomial};
/// use num::{BigRational, BigInt};
/// let v = [-2, 0, 1].iter().map(|x| BigRational::from(BigInt::from(*x))).collect::<Vec<_>>();
/// let f = Polynomial::new(v); // f=x^2-2
/// let limit = BigRational::new(BigInt::from(1), BigInt::from(10)); // 1/10
/// let left = BigRational::from(BigInt::from(0));
/// let right = BigRational::from(BigInt::from(2));
/// // (left, right] = (0, 2]
/// let (l, r) = bisection_method::<BigRational>(&f, &limit, left, right);
/// assert!(&r - &l < limit);
/// assert_eq!(l, BigRational::new(BigInt::from(11), BigInt::from(8))); // 11/8 = 1.375
/// assert_eq!(r, BigRational::new(BigInt::from(23), BigInt::from(16))); // 23/16 = 1.4375
/// ```
pub fn bisection_method<K>(f: &Polynomial<K>, limit: &K, mut left: K, mut right: K) -> (K, K)
where
    K: Sized + Clone + Ord + Zero + One + for<'x> AddAssign<&'x K> + for<'x> MulAssign<&'x K>,
    for<'x> &'x K: Add<Output = K> + Sub<Output = K> + Mul<Output = K> + Div<Output = K>,
{
    let two = K::one() + K::one();
    let mut r = f.eval(&right);
    if r.is_zero() {
        return (right.clone(), right);
    }
    while &(&right - &left) >= limit {
        let mid = &(&left + &right) / &two;
        let m = f.eval(&mid);
        if m.is_zero() {
            return (mid.clone(), mid);
        } else if &m * &r < K::zero() {
            left = mid;
        } else {
            right = mid;
            r = m;
        }
    }
    (left, right)
}

/// calculate Sturm chain and find all real roots by Sturm method.
///
/// ```
/// use find_real_roots_of_polynomial::SturmChain;
/// use polynomial_ring::{polynomial, Polynomial};
/// use num::{BigRational, BigInt};
/// let v = [-2, 0, 1].iter().map(|x| BigRational::from(BigInt::from(*x))).collect::<Vec<_>>();
/// let f = Polynomial::new(v); // f=x^2-2=(x-√2)(x+√2)
/// let limit = BigRational::new(BigInt::from(1), BigInt::from(16)); // 1/16
/// let sc = SturmChain::<BigRational>::new(f);
/// let roots = sc.find_all_real_roots(&limit);
/// assert_eq!(roots.len(), 2);
/// assert!(&roots[0].1 - &roots[0].0 < limit);
/// assert!(&roots[1].1 - &roots[1].0 < limit);
/// // -√2 ∈ [-93/64, -45/32] = [-1.453125, -1.40625]
/// //  √2 ∈ [45/32, 93/64] = [1.40625, 1.453125]
/// assert_eq!(roots[0].0, BigRational::new(BigInt::from(-93), BigInt::from(64)));
/// assert_eq!(roots[0].1, BigRational::new(BigInt::from(-45), BigInt::from(32)));
/// assert_eq!(roots[1].0, BigRational::new(BigInt::from(45), BigInt::from(32)));
/// assert_eq!(roots[1].1, BigRational::new(BigInt::from(93), BigInt::from(64)));
/// ```
pub struct SturmChain<K> {
    chain: Vec<Polynomial<K>>,
}

impl<K: Sized> SturmChain<K> {
    /// make Sturm chain
    ///
    /// Input polynomial is converted to square-free polynomial.
    /// Then Sturm chain calcurated.
    pub fn new(poly: Polynomial<K>) -> Self
    where
        K: Clone
            + Ord
            + Zero
            + One
            + Neg<Output = K>
            + for<'x> AddAssign<&'x K>
            + for<'x> SubAssign<&'x K>
            + for<'x> MulAssign<&'x K>
            + for<'x> DivAssign<&'x K>,
        for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
    {
        let poly = poly.square_free();
        let mut chain = vec![
            keep_sign_monic::<K>(poly.clone()),
            keep_sign_monic::<K>(poly.derivative()),
        ];
        loop {
            let n = chain.len();
            let g = -(&chain[n - 2] % &chain[n - 1]);
            if g.is_zero() {
                break;
            }
            chain.push(keep_sign_monic::<K>(g));
        }
        Self { chain }
    }
    fn count_sign_change(&self, x: &K) -> usize
    where
        K: Clone
            + Ord
            + Zero
            + Mul<Output = K>
            + for<'x> AddAssign<&'x K>
            + for<'x> MulAssign<&'x K>,
        for<'x> &'x K: Mul<Output = K>,
    {
        let v = self
            .chain
            .iter()
            .map(|p| p.eval(x))
            .filter(|v| !v.is_zero())
            .collect::<Vec<_>>();
        v.windows(2).filter(|e| &e[0] * &e[1] < K::zero()).count()
    }
    fn num_real_root(&self, left: &K, right: &K) -> usize
    where
        K: Clone
            + Ord
            + Zero
            + Mul<Output = K>
            + for<'x> AddAssign<&'x K>
            + for<'x> MulAssign<&'x K>,
        for<'x> &'x K: Mul<Output = K>,
    {
        let l = self.count_sign_change(left);
        let r = self.count_sign_change(right);
        l - r
    }
    /// This function returns half-open intervals. Each of them has exact one root.
    ///
    /// Let polynomial has $`n`$ real roots in `(left, right]`. Let the roots be
    /// $`r_0, r_1, \dots, r_{n-1}`$. Where, $`r_0 < r_1 < \dots < r_{n-1}`$. (Note. input
    /// ploynomial was converted to square-free polynomial. So, all root are simple root.)
    /// This function returns $`n`$ half-open intervals.  Let this intervals
    /// $`(l_0, u_0], (l_1, u_1], \dots, (l_{n-1}, u_{n-1}]`$. $`\forall i, r_i \in (l_i, u_i]`$.
    /// ```
    /// use find_real_roots_of_polynomial::SturmChain;
    /// use polynomial_ring::{polynomial, Polynomial};
    /// use num::{BigRational, BigInt, Zero, One};
    /// let v = [0, -1, 0, 1].iter().map(|x| BigRational::from(BigInt::from(*x))).collect::<Vec<_>>();
    /// let f = Polynomial::new(v); // f=x^3-1=(x-1)x(x+1)
    /// let m_two = BigRational::from(BigInt::from(-2));
    /// let p_two = BigRational::from(BigInt::from(2));
    /// let sc = SturmChain::<BigRational>::new(f);
    /// let intervals = sc.strum_method(&m_two, &p_two);
    /// assert_eq!(intervals.len(), 3);
    /// // -1 ∈ (-2, -1]
    /// assert_eq!(intervals[0].0, m_two);
    /// assert_eq!(intervals[0].1, -BigRational::one());
    /// // 0 ∈ (-1, 0]
    /// assert_eq!(intervals[1].0, -BigRational::one());
    /// assert_eq!(intervals[1].1, BigRational::zero());
    /// // 1 ∈ (0, 2]
    /// assert_eq!(intervals[2].0, BigRational::zero());
    /// assert_eq!(intervals[2].1, p_two);
    /// ```
    pub fn strum_method(&self, left: &K, right: &K) -> Vec<(K, K)>
    where
        K: Clone
            + Ord
            + Zero
            + One
            + Mul<Output = K>
            + Div<Output = K>
            + for<'x> AddAssign<&'x K>
            + for<'x> MulAssign<&'x K>,
        for<'x> &'x K: Add<Output = K> + Mul<Output = K>,
    {
        let n = self.num_real_root(left, right);
        if n == 0 {
            return Vec::new();
        } else if n == 1 {
            return vec![(left.clone(), right.clone())];
        }
        let two = K::one() + K::one();
        let mid = (left + right) / two;
        let mut v1 = self.strum_method(left, &mid);
        let mut v2 = self.strum_method(&mid, right);
        v1.append(&mut v2);
        v1
    }
    /// Find all real roots.
    ///
    /// size of interval is less than `limit`.
    /// ```
    /// use find_real_roots_of_polynomial::SturmChain;
    /// use polynomial_ring::{polynomial, Polynomial};
    /// use num::{BigRational, BigInt};
    /// let v = [-1, -1, 1].iter().map(|x| BigRational::from(BigInt::from(*x))).collect::<Vec<_>>();
    /// let f = Polynomial::new(v); // f=x^2-x-1=(x-(1+√5)/2)(x-(1-√5)/2)
    /// let limit = BigRational::new(BigInt::from(1), BigInt::from(8)); // 1/8
    /// let sc = SturmChain::<BigRational>::new(f);
    /// let roots = sc.find_all_real_roots(&limit);
    /// assert_eq!(roots.len(), 2);
    /// assert!(&roots[0].1 - &roots[0].0 < limit);
    /// assert!(&roots[1].1 - &roots[1].0 < limit);
    /// // (1-√5)/2(≒-0.61803) ∈ [-5/8, -9/16] = [-0.625, -0.5625]
    /// // (1+√5)/2(≒ 1.61803) ∈ [45/32, 3/2] = [1.5625, 1.625]
    /// assert_eq!(roots[0].0, BigRational::new(BigInt::from(-5), BigInt::from(8)));
    /// assert_eq!(roots[0].1, BigRational::new(BigInt::from(-9), BigInt::from(16)));
    /// assert_eq!(roots[1].0, BigRational::new(BigInt::from(25), BigInt::from(16)));
    /// assert_eq!(roots[1].1, BigRational::new(BigInt::from(13), BigInt::from(8)));
    /// ```
    pub fn find_all_real_roots(&self, limit: &K) -> Vec<(K, K)>
    where
        K: Clone
            + Ord
            + Zero
            + One
            + Mul<Output = K>
            + Div<Output = K>
            + Neg<Output = K>
            + for<'x> AddAssign<&'x K>
            + for<'x> MulAssign<&'x K>
            + for<'x> DivAssign<&'x K>,
        for<'x> &'x K:
            Add<Output = K> + Sub<Output = K> + Neg<Output = K> + Mul<Output = K> + Div<Output = K>,
    {
        let f = &self.chain[0];
        let m = get_range(f.clone());
        let v = self.strum_method(&-&m, &m);
        v.into_iter()
            .map(|(l, r)| bisection_method(f, limit, l, r))
            .collect::<Vec<_>>()
    }
    /// Determine if input polynomial has an imaginary root.
    ///
    /// ```
    /// use find_real_roots_of_polynomial::SturmChain;
    /// use polynomial_ring::{polynomial, Polynomial};
    /// use num::{BigRational, BigInt};
    /// let v0 = [-1, 0, 1].iter().map(|x| BigRational::from(BigInt::from(*x))).collect::<Vec<_>>();
    /// let v1 = [1, 0, 1].iter().map(|x| BigRational::from(BigInt::from(*x))).collect::<Vec<_>>();
    /// let f = Polynomial::new(v0); // f=x^2-1=(x+1)(x-1)
    /// let g = Polynomial::new(v1); // g=x^2+1=(x+i)(x-i)
    /// let f_sc = SturmChain::<BigRational>::new(f);
    /// let g_sc = SturmChain::<BigRational>::new(g);
    /// assert_eq!(f_sc.has_imaginary_root(), false);
    /// assert_eq!(g_sc.has_imaginary_root(), true);
    /// ```
    pub fn has_imaginary_root(&self) -> bool
    where
        K: Clone
            + Ord
            + Zero
            + One
            + Mul<Output = K>
            + Neg<Output = K>
            + for<'x> AddAssign<&'x K>
            + for<'x> MulAssign<&'x K>
            + for<'x> DivAssign<&'x K>,
        for<'x> &'x K: Neg<Output = K> + Mul<Output = K>,
    {
        let f = &self.chain[0];
        if let Some(d) = f.deg() {
            let m = get_range(f.clone());
            let n = self.num_real_root(&-&m, &m);
            debug_assert!(n <= d);
            n < d
        } else {
            // f = 0
            false
        }
    }
}

#[test]
fn random_test() {
    use num::{BigInt, BigRational};
    use polynomial_ring::Polynomial;
    use rand::distributions::Uniform;
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    let side = Uniform::new(-10, 10);
    for _ in 0..100 {
        let mut v = Vec::new();
        for _ in 0..10 {
            v.push(rng.sample(side));
        }
        let v = v
            .iter()
            .map(|x| BigRational::from(BigInt::from(*x)))
            .collect::<Vec<_>>();
        let f = Polynomial::new(v);
        let f = f.square_free();
        let limit = BigRational::new(BigInt::from(1), BigInt::from(8));
        let sc = SturmChain::<BigRational>::new(f.clone());
        let roots = sc.find_all_real_roots(&limit);
        for r in roots {
            assert!(&r.1 - &r.0 < limit);
            let l = f.eval(&r.0);
            let u = f.eval(&r.1);
            assert!(l * u <= BigRational::zero());
        }
    }
}

#[test]
fn a() {
    use num::{BigInt, BigRational, One};
    use polynomial_ring::Polynomial;
    let mut f = Polynomial::<BigRational>::one();
    for i in 1..=20 {
        let v = [-i, 1]
            .iter()
            .map(|x| BigRational::from(BigInt::from(*x)))
            .collect::<Vec<_>>();
        let p = Polynomial::new(v);
        f *= p;
    }
    let f = f;
    let limit = BigRational::new(BigInt::from(1), BigInt::from(8));
    let sc = SturmChain::<BigRational>::new(f);
    let roots = sc.find_all_real_roots(&limit);
    for i in 1..=20 {
        let r = BigRational::from(BigInt::from(i));
        let (l, u) = &roots[i - 1];
        assert!(l <= &r);
        assert!(u >= &r);
    }
}
