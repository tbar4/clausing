use rand::Rng;
use std::f64::consts::PI;

/// Structure to hold input parameters for the Clausing factor calculation
#[derive(Debug)]
pub struct ClausingParams {
    pub thick_screen: f64,
    pub thick_accel: f64,
    pub r_screen: f64,
    pub r_accel: f64,
    pub grid_space: f64,
    pub npart: usize,
}

/// Structure to hold the results of the Clausing calculation
#[derive(Debug)]
pub struct ClausingResults {
    pub clausing_factor: f64,
    pub max_count: usize,
    pub nlost: usize,
    pub den_cor: f64, // Downstream correction factor
}

/// Monte Carlo routine that calculates Clausing factor for CEX
/// Returns Clausing Factor and Downstream Correction factor
pub fn clausing(params: ClausingParams) -> ClausingResults {
    let mut rng = rand::rng();

    // Calculate normalized dimensions (assumes rTop = 1)
    let r_bottom = params.r_screen / params.r_accel;
    let len_bottom = (params.thick_screen + params.grid_space) / params.r_accel;
    let len_top = params.thick_accel / params.r_accel;
    let length = len_top + len_bottom;

    let mut iescape = 0;
    let mut max_count = 0;
    let mut nlost = 0;
    let mut vztot = 0.0;
    let mut vz0tot = 0.0;

    // Main particle loop
    for _ipart in 0..params.npart {
        // Launch from bottom
        let mut r0 = r_bottom * rng.random::<f64>().sqrt();
        let mut z0 = 0.0;

        let mut costheta = (1.0 - rng.random::<f64>()).sqrt();
        if costheta > 0.99999 {
            costheta = 0.99999;
        }

        let phi = 2.0 * PI * rng.random::<f64>();
        let sintheta = (1.0 - costheta.powi(2)).sqrt();

        let mut vx = phi.cos() * sintheta;
        let mut vy = phi.sin() * sintheta;
        let mut vz = costheta;

        let mut rf = r_bottom;
        let mut t = (vx * r0 + ((vx.powi(2) + vy.powi(2)) * rf.powi(2) - (vy * r0).powi(2)).sqrt())
            / (vx.powi(2) + vy.powi(2));
        let mut z = z0 + vz * t;

        vz0tot += vz;

        let mut icount = 0;
        let mut notgone = true;

        while notgone {
            icount += 1;

            // Hit wall of bottom cylinder and is re-emitted
            if z < len_bottom {
                r0 = r_bottom;
                z0 = z;

                costheta = (1.0 - rng.random::<f64>()).sqrt();
                if costheta > 0.99999 {
                    costheta = 0.99999;
                }

                let phi = 2.0 * PI * rng.random::<f64>();
                let sintheta = (1.0 - costheta.powi(2)).sqrt();

                vz = phi.cos() * sintheta;
                vy = phi.sin() * sintheta;
                vx = costheta;

                rf = r_bottom;
                t = (vx * r0 + ((vx.powi(2) + vy.powi(2)) * rf.powi(2) - (vy * r0).powi(2)).sqrt())
                    / (vx.powi(2) + vy.powi(2));
                z = z0 + t * vz;
            }

            // Emitted below but going up
            if z >= len_bottom && z0 < len_bottom {
                // Find radius at len_bottom
                t = (len_bottom - z0) / vz;
                let r = ((r0 - vx * t).powi(2) + (vy * t).powi(2)).sqrt();

                if r <= 1.0 {
                    // Continuing upward
                    rf = 1.0;
                    t = (vx * r0
                        + ((vx.powi(2) + vy.powi(2)) * rf.powi(2) - (vy * r0).powi(2)).sqrt())
                        / (vx.powi(2) + vy.powi(2));
                    z = z0 + vz * t;
                } else {
                    // Hit the upstream side of the accel grid and is re-emitted downward
                    r0 = r;
                    z0 = len_bottom;

                    costheta = (1.0 - rng.random::<f64>()).sqrt();
                    if costheta > 0.99999 {
                        costheta = 0.99999;
                    }

                    let phi = 2.0 * PI * rng.random::<f64>();
                    let sintheta = (1.0 - costheta.powi(2)).sqrt();

                    vx = phi.cos() * sintheta;
                    vy = phi.sin() * sintheta;
                    vz = -costheta;

                    rf = r_bottom;
                    t = (vx * r0
                        + ((vx.powi(2) + vy.powi(2)) * rf.powi(2) - (vy * r0).powi(2)).sqrt())
                        / (vx.powi(2) + vy.powi(2));
                    z = z0 + vz * t;
                }
            }

            // Hit the upper cylinder wall and is re-emitted
            if z >= len_bottom && z <= length {
                r0 = 1.0;
                z0 = z;

                costheta = (1.0 - rng.random::<f64>()).sqrt();
                if costheta > 0.99999 {
                    costheta = 0.99999;
                }

                let phi = 2.0 * PI * rng.random::<f64>();
                let sintheta = (1.0 - costheta.powi(2)).sqrt();

                vz = phi.cos() * sintheta;
                vy = phi.sin() * sintheta;
                vx = costheta;

                rf = 1.0;
                t = (vx * r0 + ((vx.powi(2) + vy.powi(2)) * rf.powi(2) - (vy * r0).powi(2)).sqrt())
                    / (vx.powi(2) + vy.powi(2));
                z = z0 + t * vz;

                // Find z when particle hits the bottom cylinder
                if z < len_bottom {
                    rf = r_bottom;
                    let discriminant = (vx.powi(2) + vy.powi(2)) * rf.powi(2) - (vy * r0).powi(2);

                    if discriminant < 0.0 {
                        // If sqrt argument is less than 0 then set sqrt term to 0
                        t = (vx * r0) / (vx.powi(2) + vy.powi(2));
                    } else {
                        t = (vx * r0 + discriminant.sqrt()) / (vx.powi(2) + vy.powi(2));
                    }
                    z = z0 + vz * t;
                }
            }

            // Check exit conditions
            if z < 0.0 {
                notgone = false;
            }

            if z > length {
                iescape += 1;
                vztot += vz;
                notgone = false;
            }

            if icount > 1000 {
                notgone = false;
                nlost += 1;
            }

            if icount > max_count {
                max_count = icount;
            }
        }
    }

    // Calculate results
    let clausing_factor = (r_bottom.powi(2) * iescape as f64) / params.npart as f64;
    let vz0av = vz0tot / params.npart as f64;
    let vzav = vztot / iescape as f64;
    let den_cor = vz0av / vzav; // Downstream correction factor

    ClausingResults {
        clausing_factor,
        max_count,
        nlost,
        den_cor,
    }
}

fn main() {
    // Example usage with sample parameters
    let params = ClausingParams {
        thick_screen: 1.0,
        thick_accel: 0.5,
        r_screen: 2.0,
        r_accel: 1.0,
        grid_space: 0.3,
        npart: 10000,
    };

    println!("Running Clausing factor calculation...");
    println!("Parameters: {:?}", params);
    println!();

    let results = clausing(params);

    println!("Results:");
    println!("  Clausing Factor: {:.6}", results.clausing_factor);
    println!("  Max Count: {}", results.max_count);
    println!("  Particles Lost: {}", results.nlost);
    println!("  Downstream Correction Factor: {:.6}", results.den_cor);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clausing_basic() {
        let params = ClausingParams {
            thick_screen: 1.0,
            thick_accel: 0.5,
            r_screen: 2.0,
            r_accel: 1.0,
            grid_space: 0.3,
            npart: 1000,
        };

        let results = clausing(params);

        // Basic sanity checks
        assert!(results.clausing_factor > 0.0);
        assert!(results.clausing_factor <= 4.0); // r_bottom^2 = 4.0 is theoretical max
        assert!(results.den_cor > 0.0);
    }
}
