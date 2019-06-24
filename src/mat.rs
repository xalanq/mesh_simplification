use crate::{Flt, Vct, EPS};
use std::fmt;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign};

#[derive(Copy, Clone, PartialEq, Debug, Default)]
pub struct Mat {
    pub data: [[Flt; 4]; 4],
}

impl fmt::Display for Mat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut res = String::new();
        for i in 0..4 {
            for j in 0..4 {
                res += &format!("{}{}", self[i][j], if j == 3 { "" } else { "  " });
            }
            res += if i == 3 { "" } else { "\n" };
        }
        write!(f, "{}", res)
    }
}

impl Index<usize> for Mat {
    type Output = [Flt; 4];

    fn index(&self, idx: usize) -> &[Flt; 4] {
        &self.data[idx]
    }
}

impl IndexMut<usize> for Mat {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.data[idx]
    }
}

impl Add<Mat> for Mat {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut ret = self.clone();
        ret += rhs;
        ret
    }
}

impl AddAssign<Mat> for Mat {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..4 {
            for j in 0..4 {
                self[i][j] += rhs[i][j];
            }
        }
    }
}

impl Sub<Mat> for Mat {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut ret = self.clone();
        ret -= rhs;
        ret
    }
}

impl SubAssign<Mat> for Mat {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..4 {
            for j in 0..4 {
                self[i][j] -= rhs[i][j];
            }
        }
    }
}

impl Mul<Mat> for Mat {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let mut ret = Self::default();
        for i in 0..4 {
            for j in 0..4 {
                let mut tmp = 0.0;
                for k in 0..4 {
                    tmp += self[i][k] * rhs[k][j];
                }
                ret[i][j] = tmp;
            }
        }
        ret
    }
}

impl Mat {
    pub fn identity() -> Self {
        return Self {
            data: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };
    }

    pub fn multiply_by_vct(&self, v: Vct) -> Flt {
        // vQv^T
        let mut tmp = [0.0, 0.0, 0.0, 0.0];
        for i in 0..4 {
            for j in 0..4 {
                tmp[i] += v[j] * self[j][i];
            }
        }
        let mut ret = 0.0;
        for i in 0..4 {
            ret += tmp[i] * v[i];
        }
        ret
    }

    pub fn split(&self) -> Self {
        let mut ret = self.clone();
        ret[3][0] = 0.0;
        ret[3][1] = 0.0;
        ret[3][2] = 0.0;
        ret[3][3] = 1.0;
        ret
    }

    pub fn inverse(&self) -> Option<Self> {
        let mut a = self.clone();
        let mut b = Self::identity();
        for x in 0..4 {
            let mut z = x;
            for i in x + 1..4 {
                if a[i][x].abs() > a[z][x].abs() {
                    z = i;
                }
            }
            if a[z][x].abs() <= EPS {
                return None;
            }
            if z != x {
                for y in x..4 {
                    let t = a[z][y];
                    a[z][y] = a[x][y];
                    a[x][y] = t;
                }
                for y in 0..4 {
                    let t = b[z][y];
                    b[z][y] = b[x][y];
                    b[x][y] = t;
                }
            }
            let inv = -1.0 / a[x][x];
            for i in 0..4 {
                if i == x {
                    continue;
                }
                let d = inv * a[i][x];
                for y in x..4 {
                    a[i][y] += d * a[x][y];
                }
                for y in 0..4 {
                    b[i][y] += d * b[x][y];
                }
            }
        }
        for x in 0..4 {
            let inv = 1.0 / a[x][x];
            for y in 0..4 {
                b[x][y] *= inv;
            }
        }
        Some(b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse() {
        let a = Mat {
            data: [
                /*
                [0.326, 0.274, 0.526, 2.267],
                [1.916, 0.195, 0.839, 0.649],
                [0.458, 1.867, 0.476, 0.987],
                [0.996, 0.507, 1.616, 0.927],
                */
                [1.906, -4.297, -0.778, 9.577],
                [-4.297, 9.755, 1.743, -21.651],
                [-0.778, 1.743, 0.339, -3.865],
                [0.000, 0.000, 0.000, 1.000],
            ],
        };
        let b = a.inverse().unwrap();
        let c = a * b;
        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    assert!((c[i][j] - 1.0).abs() < EPS);
                } else {
                    assert!(c[i][j].abs() < EPS);
                }
            }
        }
    }

}
