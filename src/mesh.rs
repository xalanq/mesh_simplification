use crate::{Flt, Mat, Vct};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

const COST_EPS: Flt = 1e50;
const DIST_EPS: Flt = 1e50;

#[derive(Clone, Copy, Debug, PartialEq)]
struct Edge {
    pub v1: usize,
    pub v2: usize,
    pub qv: Mat,
    pub v: Vct,
    pub cost: Flt,
}

impl Edge {
    pub fn new(v1: usize, v2: usize, qv: Mat, v: Vct, cost: Flt) -> Self {
        Self { v1, v2, qv, v, cost }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct State {
    pub cost: Flt,
    pub index: usize,
}

impl State {
    pub fn new(cost: Flt, index: usize) -> Self {
        Self { cost, index }
    }
}

impl Eq for State {}

impl Ord for State {
    fn cmp(&self, other: &State) -> Ordering {
        other.partial_cmp(self).unwrap()
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &State) -> Option<Ordering> {
        other.cost.partial_cmp(&self.cost)
    }
}

type Tri = (usize, usize, usize);

#[derive(Clone, Debug)]
pub struct Mesh {
    pub pos: Vec<Vct>,
    pub tri: Vec<Tri>,
}

impl Mesh {
    pub fn new(path: &str) -> Self {
        let (pos, tri) = Self::load(path);
        Self { pos, tri }
    }

    fn cal_q(v1: &Vct, v2: &Vct, v3: &Vct) -> Mat {
        let mut q = Mat::default();
        let norm = ((*v1 - *v3) % (*v2 - *v3)).norm();
        let v = [norm.x, norm.y, norm.z, -norm.dot(*v3)];
        for i in 0..4 {
            for j in 0..4 {
                q[i][j] = v[i] * v[j];
            }
        }
        q
    }

    pub fn simplify(&self, ratio: Flt) -> Self {
        let (mut pos, mut tri) = (self.pos.clone(), self.tri.clone());
        let mut q = vec![Mat::default(); pos.len()];
        let mut qp = Vec::with_capacity(tri.len());
        let mut head = vec![vec![0; 0]; pos.len()];

        tri.iter().enumerate().for_each(|(i, &(v1, v2, v3))| {
            let qk = Self::cal_q(&pos[v1], &pos[v2], &pos[v3]);
            qp.push(qk);
            q[v1] += qk;
            q[v2] += qk;
            q[v3] += qk;
            head[v1].push(i);
            head[v2].push(i);
            head[v3].push(i);
        });

        let mut dele = vec![false; pos.len()];
        let mut edge = vec![];
        let mut heap = BinaryHeap::new();

        macro_rules! test_edge {
            ($v1:expr, $v2:expr) => {
                let e = pos[$v1] - pos[$v2];
                if e.len2() >= DIST_EPS {
                    return;
                }
                let qv = q[$v1] + q[$v2];
                let inv = qv.split().inverse();
                let v = match inv {
                    Some(inv) => Vct::new(inv[0][3], inv[1][3], inv[2][3]),
                    None => (pos[$v1] + pos[$v2]) * 0.5,
                };
                let cost = qv.multiply_by_vct(v);
                if cost >= COST_EPS {
                    return;
                }
                edge.push(Edge::new($v1, $v2, qv, v, cost));
                heap.push(State::new(cost, edge.len() - 1));
            };
        }

        macro_rules! in_tri {
            ($i:expr, $v:expr) => {
                tri[$i].0 == $v || tri[$i].1 == $v || tri[$i].2 == $v
            };
        }

        macro_rules! is_valid_tri {
            ($i:expr) => {
                !(dele[tri[$i].0] || dele[tri[$i].1] || dele[tri[$i].2])
            };
        }

        tri.iter().for_each(|&(v1, v2, v3)| {
            test_edge!(v1, v2);
            test_edge!(v2, v3);
            test_edge!(v1, v3);
        });

        let mut limit = (self.tri.len() as Flt * (1.0 - ratio)) as i64;
        // 一次缩掉两个面
        while let Some(state) = heap.pop() {
            let e = edge[state.index];
            if dele[e.v1] || dele[e.v2] {
                continue;
            }
            let mut head_v = vec![];
            let mut edge_v = vec![];
            q.push(Mat::default());
            pos.push(e.v);
            dele.push(false);
            let v = pos.len() - 1;

            // 只有一个顶点是v1的三角形
            macro_rules! type_a {
                ($i:expr, $v1:expr) => {
                    if tri[$i].1 == $v1 {
                        tri[$i] = (tri[$i].1, tri[$i].2, tri[$i].0);
                    } else if tri[$i].2 == $v1 {
                        tri[$i] = (tri[$i].2, tri[$i].0, tri[$i].1);
                    }
                    let (v2, v3) = (tri[$i].1, tri[$i].2);
                    let qk = Self::cal_q(&pos[v], &pos[v2], &pos[v3]);
                    let dq = qk - qp[$i];
                    head_v.push($i);
                    q[v] += qk;
                    q[v2] += dq;
                    q[v3] += dq;
                    qp[$i] = qk;
                    tri[$i].0 = v;
                    edge_v.push(v2);
                    edge_v.push(v3);
                };
            }

            // 有两个顶点是v1和v2的三角形
            macro_rules! type_b {
                ($i:expr, $v1:expr, $v2:expr) => {
                    let v3 = tri[$i].0 + tri[$i].1 + tri[$i].2 - $v1 - $v2;
                    q[v3] -= qp[$i];
                    edge_v.push(v3);
                };
            }

            head[e.v1].iter().for_each(|&i| {
                if is_valid_tri!(i) {
                    if !in_tri!(i, e.v2) {
                        type_a!(i, e.v1);
                    } else {
                        type_b!(i, e.v1, e.v2);
                    }
                }
            });
            head[e.v2].iter().for_each(|&i| {
                if is_valid_tri!(i) && !in_tri!(i, e.v1) {
                    type_a!(i, e.v2);
                }
            });
            dele[e.v1] = true;
            dele[e.v2] = true;
            edge_v.iter().for_each(|&vi| {
                if !dele[vi] {
                    test_edge!(v, vi);
                    dele[vi] = true;
                }
            });
            edge_v.iter().for_each(|&vi| {
                dele[vi] = false;
            });
            head.push(head_v);
            limit -= 2;
            if limit < 0 {
                break;
            }
        }

        let (mut new_pos, mut new_tri) = (vec![], vec![]);
        let mut cnt = 0;
        let mut id = vec![-1 as i64; pos.len()];
        macro_rules! gg {
            ($i:expr) => {{
                if id[$i] == -1 {
                    id[$i] = cnt;
                    new_pos.push(pos[$i]);
                    cnt += 1;
                }
                id[$i] as usize
            }};
        }
        for i in 0..tri.len() {
            if is_valid_tri!(i) {
                new_tri.push((gg!(tri[i].0), gg!(tri[i].1), gg!(tri[i].2)));
            }
        }
        Self { pos: new_pos, tri: new_tri }
    }

    fn load(path: &str) -> (Vec<Vct>, Vec<(usize, usize, usize)>) {
        println!("Loading the object from {}", path);
        let file = File::open(path).expect(&format!("Cannot open {}", path));
        let (mut t_v, mut t_f) = (vec![], vec![]);
        for line in BufReader::new(file).lines() {
            let line = line.expect("Failed to load the mesh object");
            let mut w = line.split_whitespace();
            macro_rules! nx {
                () => {
                    w.next().unwrap().parse().unwrap()
                };
            }
            macro_rules! nxtf {
                () => {{
                    let mut a = Vec::new();
                    w.next().unwrap().split('/').for_each(|x| {
                        if let Ok(i) = x.parse::<usize>() {
                            a.push(i);
                        }
                    });
                    match a.len() {
                        1 => (a[0], 0, 0),
                        2 => (a[0], 0, a[1]),
                        3 => (a[0], a[1], a[2]),
                        _ => panic!("invalid vertex of a face"),
                    }
                }};
            }
            macro_rules! wp {
                ($e:expr) => {{
                    $e;
                    w.next().map(|_| panic!("The mesh object has a non-triangle"));
                }};
            }
            match w.next() {
                Some("v") => wp!(t_v.push(Vct::new(nx!(), nx!(), nx!()))),
                Some("f") => wp!(t_f.push((nxtf!(), nxtf!(), nxtf!()))),
                _ => (),
            }
        }
        let mut vis = HashMap::new();
        let (mut pos, mut tri) = (vec![], vec![]);
        macro_rules! gg {
            ($a:expr) => {{
                *vis.entry($a).or_insert_with(|| {
                    pos.push(t_v[$a.0 - 1]);
                    pos.len() - 1
                })
            }};
        }
        t_f.iter().for_each(|&(a, b, c)| {
            let g = (gg!(a), gg!(b), gg!(c));
            tri.push(g);
        });
        println!("...Loaded");
        (pos, tri)
    }

    pub fn save(&self, path: &str) {
        println!("Saving the object to {}", path);
        let mut file = File::create(path).expect(&format!("Cannot open {}", path));
        let mut s = String::new();
        self.pos.iter().for_each(|p| {
            s += &format!("v {} {} {}\n", p.x, p.y, p.z);
        });
        self.tri.iter().for_each(|&(a, b, c)| {
            s += &format!("f {} {} {}\n", a + 1, b + 1, c + 1);
        });
        write!(file, "{}", s).expect(&format!("Cannot write to file"));
        println!("...Saved");
    }
}
