calculate Sturm chain and find all real roots by Sturm method.

```rust
use find_real_roots_of_polynomial::SturmChain;
use polynomial_ring::{polynomial, Polynomial};
use num::{BigRational, BigInt};
let v = [-2, 0, 1].iter().map(|x| BigRational::from(BigInt::from(*x))).collect::<Vec<_>>();
let f = Polynomial::new(v); // f=x^2-2=(x-√2)(x+√2)
let limit = BigRational::new(BigInt::from(1), BigInt::from(16)); // 1/16
let sc = SturmChain::<BigRational>::new(f);
let roots = sc.find_all_real_roots(&limit);
assert_eq!(roots.len(), 2);
assert!(&roots[0].1 - &roots[0].0 < limit);
assert!(&roots[1].1 - &roots[1].0 < limit);
// -√2 ∈ [-93/64, -45/32] = [-1.453125, -1.40625]
//  √2 ∈ [45/32, 93/64] = [1.40625, 1.453125]
assert_eq!(roots[0].0, BigRational::new(BigInt::from(-93), BigInt::from(64)));
assert_eq!(roots[0].1, BigRational::new(BigInt::from(-45), BigInt::from(32)));
assert_eq!(roots[1].0, BigRational::new(BigInt::from(45), BigInt::from(32)));
assert_eq!(roots[1].1, BigRational::new(BigInt::from(93), BigInt::from(64)));
```
