use crate::arith::U256;
use crate::fields::{const_fq, FieldElement, Fq, Fq2, Fq6};
use core::ops::{Add, Div, Mul, Neg, Sub};
use rand::Rng;

fn frobenius_coeffs_c1(power: usize) -> Fq2 {
    match power % 12 {
        0 => Fq2::one(),
        1 => Fq2::new(
            const_fq([
                0xd60b35dadcc9e470,
                0x5c521e08292f2176,
                0xe8b99fdd76e68b60,
                0x1284b71c2865a7df,
            ]),
            const_fq([
                0xca5cf05f80f362ac,
                0x747992778eeec7e5,
                0xa6327cfe12150b8e,
                0x246996f3b4fae7e6,
            ]),
        ),
        2 => Fq2::new(
            const_fq([
                0xe4bd44e5607cfd49,
                0xc28f069fbb966e3d,
                0x5e6dd9e7e0acccb0,
                0x30644e72e131a029,
            ]),
            Fq::zero(),
        ),
        3 => Fq2::new(
            const_fq([
                0xe86f7d391ed4a67f,
                0x894cb38dbe55d24a,
                0xefe9608cd0acaa90,
                0x19dc81cfcc82e4bb,
            ]),
            const_fq([
                0x7694aa2bf4c0c101,
                0x7f03a5e397d439ec,
                0x6cbeee33576139d,
                0xabf8b60be77d73,
            ]),
        ),
        _ => unimplemented!(),
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq12 {
    c0: Fq6,
    c1: Fq6,
}

impl Fq12 {
    pub fn new(c0: Fq6, c1: Fq6) -> Self {
        Fq12 { c0, c1 }
    }

    fn final_exponentiation_first_chunk(&self) -> Option<Fq12> {
        match self.inverse() {
            Some(b) => {
                let a = self.unitary_inverse();
                let c = a * b;
                let d = c.frobenius_map(2);

                Some(d * c)
            }
            None => None,
        }
    }

    fn final_exponentiation_last_chunk(&self) -> Fq12 {
        let a = self.exp_by_neg_z();
        let b = a.cyclotomic_squared();
        let c = b.cyclotomic_squared();
        let d = c * b;

        let e = d.exp_by_neg_z();
        let f = e.cyclotomic_squared();
        let g = f.exp_by_neg_z();
        let h = d.unitary_inverse();
        let i = g.unitary_inverse();

        let j = i * e;
        let k = j * h;
        let l = k * b;
        let m = k * e;
        let n = *self * m;

        let o = l.frobenius_map(1);
        let p = o * n;

        let q = k.frobenius_map(2);
        let r = q * p;

        let s = self.unitary_inverse();
        let t = s * l;
        let u = t.frobenius_map(3);

        u * r
    }

    pub fn final_exponentiation(&self) -> Option<Fq12> {
        self.final_exponentiation_first_chunk()
            .map(|a| a.final_exponentiation_last_chunk())
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        Fq12 {
            c0: self.c0.frobenius_map(power),
            c1: self
                .c1
                .frobenius_map(power)
                .scale(frobenius_coeffs_c1(power)),
        }
    }

    pub fn exp_by_neg_z(&self) -> Fq12 {
        self.cyclotomic_pow(U256::from([4965661367192848881, 0, 0, 0]))
            .unitary_inverse()
    }

    pub fn unitary_inverse(&self) -> Fq12 {
        Fq12::new(self.c0, -self.c1)
    }

    pub fn mul_by_024(&self, ell_0: Fq2, ell_vw: Fq2, ell_vv: Fq2) -> Fq12 {
        let mut z0 = self.c0.c0;
        let mut z1 = self.c0.c1;
        let mut z2 = self.c0.c2;
        let mut z3 = self.c1.c0;
        let mut z4 = self.c1.c1;
        let mut z5 = self.c1.c2;
        let x0 = ell_0;
        let x2 = ell_vv;
        let x4 = ell_vw;

        let mut d0 = z0;
        d0.mul_inp(&x0);
        let mut d2 = z2;
        d2.mul_inp(&x2);
        let mut d4 = z4;
        d4.mul_inp(&x4);

        let mut t2 = z0;
        t2.add_inp(&z4);
        let mut t1 = z0;
        t1.add_inp(&z2);
        let mut s0 = z1;
        s0.add_inp(&z3);
        s0.add_inp(&z5);

        let mut s1 = z1;
        s1.mul_inp(&x2);
        let mut t3 = s1;
        t3.add_inp(&d4);
        let mut t4 = t3;
        t4.mul_by_nonresidue_inp();
        t4.add_inp(&d0);
        z0 = t4;

        let mut t3 = z5;
        t3.mul_inp(&x4);
        s1.add_inp(&t3);
        t3.add_inp(&d2);
        t4 = t3;
        t4.mul_by_nonresidue_inp();
        t3 = z1;
        t3.mul_inp(&x0);
        s1.add_inp(&t3);
        t4.add_inp(&t3);
        z1 = t4;

        let mut t0 = x0;
        t0.add_inp(&x2);
        t3 = t1;
        t3.mul_inp(&t0);
        t3.sub_inp(&d0);
        t3.sub_inp(&d2);
        t4 = z3;
        t4.mul_inp(&x4);
        s1.add_inp(&t4);
        t3.add_inp(&t4);
        t0 = z2;
        t0.add_inp(&z4);
        z2 = t3;

        let mut t1 = x2;
        t1.add_inp(&x4);
        t3 = t0;
        t3.mul_inp(&t1);
        t3.sub_inp(&d2);
        t3.sub_inp(&d4);
        t4 = t3;
        t4.mul_by_nonresidue_inp();
        t3 = z3;
        t3.mul_inp(&x0);
        s1.add_inp(&t3);
        t4.add_inp(&t3);
        z3 = t4;

        t3 = z5;
        t3.mul_inp(&x2);
        s1.add_inp(&t3);
        t4 = t3;
        t4.mul_by_nonresidue_inp();
        t0 = x0;
        t0.add_inp(&x4);
        t3 = t2;
        t3.mul_inp(&t0);
        t3.sub_inp(&d0);
        t3.sub_inp(&d4);
        t4.add_inp(&t3);
        z4 = t4;

        t0 = x0;
        t0.add_inp(&x2);
        t0.add_inp(&x4);
        t3 = s0;
        t3.mul_inp(&t0);
        t3.sub_inp(&s1);
        z5 = t3;

        Fq12 {
            c0: Fq6::new(z0, z1, z2),
            c1: Fq6::new(z3, z4, z5),
        }
    }

    pub fn cyclotomic_squared(&self) -> Self {
        let z0 = self.c0.c0;
        let z4 = self.c0.c1;
        let z3 = self.c0.c2;
        let z2 = self.c1.c0;
        let z1 = self.c1.c1;
        let z5 = self.c1.c2;

        let mut tmp = z0;
        tmp.mul_inp(&z1);
        let mut t0 = z0;
        t0.add_inp(&z1);
        let mut t1 = z1;
        t1.mul_by_nonresidue_inp();
        t1.add_inp(&z0);
        t0.mul_inp(&t1);
        t0.sub_inp(&tmp);
        t0.sub_inp(&tmp.mul_by_nonresidue());
        let mut t1 = tmp;
        t1.double_inp();

        let mut tmp = z2;
        tmp.mul_inp(&z3);
        let mut t2 = z2;
        t2.add_inp(&z3);
        let mut t3 = z3;
        t3.mul_by_nonresidue_inp();
        t3.add_inp(&z2);
        t2.mul_inp(&t3);
        t2.sub_inp(&tmp);
        t2.sub_inp(&tmp.mul_by_nonresidue());
        let mut t3 = tmp;
        t3.double_inp();

        let mut tmp = z4;
        tmp.mul_inp(&z5);
        let mut t4 = z4;
        t4.add_inp(&z5);
        let mut t5 = z5;
        t5.mul_by_nonresidue_inp();
        t5.add_inp(&z4);
        t4.mul_inp(&t5);
        t4.sub_inp(&tmp);
        t4.sub_inp(&tmp.mul_by_nonresidue());
        let mut t5 = tmp;
        t5.double_inp();

        let mut new_z0 = t0;
        new_z0.sub_inp(&z0);
        new_z0.double_inp();
        new_z0.add_inp(&t0);

        let mut new_z1 = t1;
        new_z1.add_inp(&z1);
        new_z1.double_inp();
        new_z1.add_inp(&t1);

        let mut new_z2 = t5;
        new_z2.mul_by_nonresidue_inp();
        new_z2.add_inp(&z2);
        new_z2.double_inp();
        new_z2.add_inp(&t5.mul_by_nonresidue());

        let mut new_z3 = t4;
        new_z3.sub_inp(&z3);
        new_z3.double_inp();
        new_z3.add_inp(&t4);

        let mut new_z4 = t2;
        new_z4.sub_inp(&z4);
        new_z4.double_inp();
        new_z4.add_inp(&t2);

        let mut new_z5 = t3;
        new_z5.add_inp(&z5);
        new_z5.double_inp();
        new_z5.add_inp(&t3);

        Fq12 {
            c0: Fq6::new(new_z0, new_z4, new_z3),
            c1: Fq6::new(new_z2, new_z1, new_z5),
        }
    }

    pub fn cyclotomic_pow<I: Into<U256>>(&self, by: I) -> Self {
        let mut res = Self::one();

        let mut found_one = false;

        for i in by.into().bits() {
            if found_one {
                res = res.cyclotomic_squared();
            }

            if i {
                found_one = true;
                res = *self * res;
            }
        }

        res
    }
}

impl FieldElement for Fq12 {
    fn zero() -> Self {
        Fq12 {
            c0: Fq6::zero(),
            c1: Fq6::zero(),
        }
    }

    fn one() -> Self {
        Fq12 {
            c0: Fq6::one(),
            c1: Fq6::zero(),
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq12 {
            c0: Fq6::random(rng),
            c1: Fq6::random(rng),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn squared(&self) -> Self {
        let ab = self.c0 * self.c1;

        Fq12 {
            c0: (self.c1.mul_by_nonresidue() + self.c0) * (self.c0 + self.c1)
                - ab
                - ab.mul_by_nonresidue(),
            c1: ab + ab,
        }
    }

    fn inverse(self) -> Option<Self> {
        (self.c0.squared() - (self.c1.squared().mul_by_nonresidue()))
            .inverse()
            .map(|t| Fq12 {
                c0: self.c0 * t,
                c1: -(self.c1 * t),
            })
    }
}

impl Mul for Fq12 {
    type Output = Fq12;

    fn mul(self, other: Fq12) -> Fq12 {
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;

        Fq12 {
            c0: bb.mul_by_nonresidue() + aa,
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - aa - bb,
        }
    }
}

impl Div for Fq12 {
    type Output = Fq12;

    fn div(self, other: Fq12) -> Fq12 {
        self * other.inverse().expect("division by zero")
    }
}

impl Sub for Fq12 {
    type Output = Fq12;

    fn sub(self, other: Fq12) -> Fq12 {
        Fq12 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
        }
    }
}

impl Add for Fq12 {
    type Output = Fq12;

    fn add(self, other: Fq12) -> Fq12 {
        Fq12 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

impl Neg for Fq12 {
    type Output = Fq12;

    fn neg(self) -> Fq12 {
        Fq12 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}
