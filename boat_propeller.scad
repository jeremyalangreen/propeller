// Parametric Boat Propeller
// ---------------------------------------------------------------
// Blades are true helical surfaces: at every radius r the blade
// follows a helix that advances `pitch` mm per revolution, so the
// local blade angle is atan(pitch / (2*PI*r)).
//
// Two blade definitions are available (see `blade_method`):
//   "simple"     - free-form outline controlled by the Simple Blade
//                  Shape parameters below.
//   "wageningen" - Wageningen B-series outline: chord, skew (leading
//                  edge position), max-thickness position, thickness,
//                  15 deg rake, hub size and B4 root pitch reduction
//                  per Kuiper (1992). Sections use a flat face with a
//                  parabolic back peaking at the B-series max-thickness
//                  position; the small face camber the series has near
//                  the edges at inner radii is not modelled.
//
// Orientation: +Z points forward (toward the boat / shaft entry).
// Thrust is produced toward +Z. Default is a right-hand prop
// (turns clockwise when viewed from astern in forward gear).
// All dimensions in mm.
// ---------------------------------------------------------------

/* [Main Dimensions] */
// Outside diameter of the prop (tip to tip)
prop_diameter = 66.675;
// Pitch: axial advance per revolution (nominal pitch at 0.7R)
pitch = 30;
// Number of blades
num_blades = 2; // [2:1:7]
// Shaft diameter (large end if the bore is tapered)
shaft_diameter = 9.525;
// Extra clearance added to the bore diameter (for 3D printing fit)
shaft_clearance = 0.2;
// Rotation direction
hand = "left"; // [right, left]

/* [Blade Method] */
// Blade geometry definition
blade_method = "simple"; // [simple, wageningen]

/* [Wageningen B-Series] */
// Expanded blade area ratio (total blade area / disc area)
expanded_area_ratio = 0.55; // [0.3:0.01:1.05]

/* [Simple Blade Shape] */
// Maximum blade chord as a fraction of the diameter
blade_width_ratio = 0.32; // [0.1:0.01:0.6]
// Blade thickness at the root as a fraction of the diameter
root_thickness_ratio = 0.045; // [0.01:0.005:0.1]
// Blade thickness at the tip (mm)
tip_thickness = 1.0;
// Skew: how far the tip sweeps back against rotation (degrees)
skew_angle = 15; // [0:1:45]
// Rake: tilt of the blades aft (degrees)
rake_angle = 0; // [-15:1:20]

/* [Blade Edges] */
// Minimum leading/trailing edge thickness (mm), applies to both methods
edge_thickness = 0.6;

/* [Hub] */
// Hub diameter at the forward end (0 = automatic)
hub_diameter = 0;
// Hub length (0 = automatic, sized to cover blade roots)
hub_length = 0;
// Aft hub diameter as a fraction of forward diameter (taper)
hub_taper = 0.85; // [0.5:0.01:1]
// Edge chamfer on the hub ends
hub_chamfer = 1;

/* [Bore] */
// Bore taper as 1:N on diameter (0 = straight bore; SAE shafts are 1:16)
bore_taper = 0;
// Keyway width (0 = no keyway)
keyway_width = 0;
// Keyway depth into hub, measured from the bore wall
keyway_depth = 0;
// Shear Pin Diameter (0 = no shear pin)
shearpin_diameter = 3.175;
// Shear Pin Length
shearpin_length = 20;

/* [Resolution] */
// Radial sections per blade
radial_steps = 32; // [8:1:80]
// Points along each chord (per side)
chord_steps = 24; // [6:1:60]
$fn = 96;

// ---------------------------------------------------------------
// Wageningen B-series data (Kuiper 1992), sampled at
// r/R = 0.2, 0.3, ... 0.9 (thickness to 1.0).
// The B3 series has its own tables; B4-B7 share one set, which is
// also used here for 2 blades.
// ---------------------------------------------------------------
B3 = num_blades == 3;

// Chord: c * Z / (D * EAR)
wb_chord = B3 ? [1.633, 1.832, 2.000, 2.120, 2.186, 2.168, 2.127, 1.657]
              : [1.662, 1.882, 2.050, 2.152, 2.187, 2.144, 1.970, 1.582];
wb_chord_max = B3 ? 2.186 : 2.187;

// Chord divided by sqrt(1 - r/R). This varies smoothly and stays finite
// at the tip, so interpolating it and multiplying back gives a rounded
// tip that closes to zero with no kink.
wb_chord_g = [ for (i = [0 : len(wb_chord) - 1])
    wb_chord[i] / sqrt(1 - (0.2 + 0.1 * i)) ];

// Leading edge distance from the generator line, as a fraction of chord
wb_lead = B3 ? [0.616, 0.611, 0.599, 0.583, 0.558, 0.526, 0.481, 0.400]
             : [0.617, 0.613, 0.601, 0.586, 0.561, 0.524, 0.463, 0.351];

// Max thickness position from the leading edge, as a fraction of chord
wb_tpos = B3 ? [0.350, 0.350, 0.350, 0.355, 0.389, 0.442, 0.478, 0.500]
             : [0.350, 0.350, 0.351, 0.355, 0.389, 0.443, 0.479, 0.500];

// Catmull-Rom interpolation of values v sampled at r/R = 0.2, 0.3, ...
// with linear extrapolation outside the table
function wb_interp(v, rR) =
    let(n = len(v), t = (rR - 0.2) / 0.1)
    t <= 0     ? v[0] + t * (v[1] - v[0])
  : t >= n - 1 ? v[n - 1] + (t - (n - 1)) * (v[n - 1] - v[n - 2])
  : let(i  = floor(t), f = t - i,
        p0 = i == 0 ? 2 * v[0] - v[1] : v[i - 1],
        p1 = v[i],
        p2 = v[i + 1],
        p3 = i + 2 >= n ? 2 * v[n - 1] - v[n - 2] : v[i + 2])
    0.5 * (2 * p1 + (p2 - p0) * f
           + (2 * p0 - 5 * p1 + 4 * p2 - p3) * f * f
           + (3 * p1 - p0 - 3 * p2 + p3) * f * f * f);

// Max thickness t / D. B4-B7: t/D = A - B * Z
wb_rR  = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
wb_A   = [0.0526, 0.0464, 0.0402, 0.0340, 0.0278, 0.0216, 0.0154, 0.0092, 0.0030];
wb_B   = [0.0040, 0.0035, 0.0030, 0.0025, 0.0020, 0.0015, 0.0010, 0.0005, 0.0000];
wb_t3  = [0.0406, 0.0359, 0.0312, 0.0265, 0.0218, 0.0171, 0.0124, 0.0077, 0.0030];
wb_thick = [ for (i = [0 : len(wb_rR) - 1])
    [wb_rR[i], B3 ? wb_t3[i] : wb_A[i] - wb_B[i] * num_blades] ];

// Pitch reduction toward the root (B4 only; other counts are constant pitch)
wb_pitch = num_blades == 4
    ? [[0.2, 0.822], [0.3, 0.887], [0.4, 0.950], [0.5, 0.992], [0.6, 1.0]]
    : [[0.2, 1.0], [1.0, 1.0]];

wb_rake = 15;

// Typical B-series ranges of EAR by blade count
wb_ear_range = num_blades == 3 ? [0.35, 0.80]
             : num_blades == 4 ? [0.40, 1.00]
             : num_blades == 5 ? [0.45, 1.05]
             : num_blades == 6 ? [0.50, 1.00]
             : num_blades == 7 ? [0.55, 0.85]
             : [0.30, 0.38];

// ---------------------------------------------------------------
// Derived values
// ---------------------------------------------------------------
WB      = blade_method == "wageningen";
R       = prop_diameter / 2;
hub_auto = WB ? prop_diameter * (B3 ? 0.18 : 0.167) : prop_diameter * 0.2;
hub_d   = hub_diameter > 0 ? hub_diameter : max(shaft_diameter * 2.2, hub_auto);
hub_r   = hub_d / 2;
hub_r_aft = hub_r * hub_taper;
H       = pitch / (2 * PI);            // nominal axial advance per radian
max_chord = WB ? wb_chord_max * prop_diameter * expanded_area_ratio / num_blades
               : blade_width_ratio * prop_diameter;
root_t  = root_thickness_ratio * prop_diameter;
r_root  = hub_r_aft * 0.75;            // blade starts buried inside the hub
rake    = WB ? wb_rake : rake_angle;

bore_d  = shaft_diameter + shaft_clearance;

if (WB) {
    if (expanded_area_ratio < wb_ear_range[0] || expanded_area_ratio > wb_ear_range[1])
        echo(str("WARNING: EAR ", expanded_area_ratio, " is outside the usual B", num_blades,
                 " range ", wb_ear_range[0], "-", wb_ear_range[1]));
    if (pitch / prop_diameter < 0.5 || pitch / prop_diameter > 1.4)
        echo(str("WARNING: P/D ", pitch / prop_diameter,
                 " is outside the B-series range 0.5-1.4"));
}

// Normalised blade position: 0 at hub surface, 1 at tip
function xr(r) = max(0, (r - hub_r) / (R - hub_r));

// --- Simple method distributions (functions of xr) ---

// Chord: moderate at the root, widest around 60%, rounded tip
function simple_chord(x) =
    max_chord * max(0.02,
        (0.55 + 0.45 * sin(min(x, 0.6) / 0.6 * 90)) * sqrt(max(0, 1 - pow(x, 3))));

// Max thickness, linear from root to tip
function simple_thick(x) = root_t + (tip_thickness - root_t) * x;

// --- Wageningen distributions (functions of r/R) ---

// Matches the table at each sampled radius and closes to zero at 1.0R
function wb_chord_at(rR) =
    prop_diameter * expanded_area_ratio / num_blades
        * wb_interp(wb_chord_g, rR) * sqrt(max(0, 1 - rR));

// --- Per-radius blade properties for the selected method ---

function chord_at(r) = max(0.02 * max_chord,
    WB ? wb_chord_at(r / R) : simple_chord(xr(r)));

function thick_at(r) = max(edge_thickness,
    WB ? prop_diameter * lookup(r / R, wb_thick) : simple_thick(xr(r)));

// Max thickness position from the leading edge (fraction of chord)
function tpos_at(r) = WB ? wb_interp(wb_tpos, r / R) : 0.5;

// Local pitch as a fraction of nominal
function pitch_fac(r) = WB ? lookup(r / R, wb_pitch) : 1;

// Angular offset (deg) of the section centre from the generator line.
// Negative = skewed back, against rotation.
function centre_off(r, c, L) =
    WB ? (wb_interp(wb_lead, r / R) - 0.5) * c / L * 180 / PI
       : -skew_angle * pow(xr(r), 1.5);

// Axial offset from rake (positive rake moves tips aft, -Z)
function rake_off(r) = -max(0, r - hub_r) * tan(rake);

// Section thickness, blended down to the edge thickness as the chord
// narrows at the tip so the section never gets thicker than it is wide
function sect_thick(r, c) =
    edge_thickness + (thick_at(r) - edge_thickness) * min(1, c / (0.3 * max_chord));

// Back-surface thickness fraction at chord position s (0 = LE, 1 = TE),
// peaking at b. Parabolic on each side of the peak.
function back_frac(s, b) =
    s < b ? 1 - pow((b - s) / b, 2) : 1 - pow((s - b) / (1 - b), 2);

// Auto hub length: cover the axial extent of the blade root plus margin
L_hub_root = sqrt(hub_r * hub_r + H * H);
root_axial = chord_at(hub_r) * H / L_hub_root + thick_at(hub_r) * hub_r / L_hub_root;
hub_len = hub_length > 0 ? hub_length
                         : max(shaft_diameter * 1.2, root_axial * 1.25 + 2 * hub_chamfer);

assert(R > hub_r + 1, "Prop diameter must be larger than the hub diameter");
assert(hub_r_aft > bore_d / 2 + 1, "Hub is too small for the shaft; increase hub_diameter");

// ---------------------------------------------------------------
// Blade geometry
// ---------------------------------------------------------------
// Point on the helical surface (advance h per radian) at radius r and
// angle phi (deg), offset n along the section normal (toward the
// forward/suction side).
function blade_pt(r, phi, n, h, L, z0) =
    let(nt = -h / L, nz = r / L)
    [ r * cos(phi) - n * nt * sin(phi),
      r * sin(phi) + n * nt * cos(phi),
      h * phi * PI / 180 + n * nz + z0 ];

K = chord_steps;
M = 2 * K;              // points per section loop
N = radial_steps;

// Radii clustered toward the tip where the outline curves fastest
function radius_at(i) = r_root + (R - r_root) * sin(i / N * 90);

function blade_points() = [
    for (i = [0 : N])
        let(r    = radius_at(i),
            c    = chord_at(r),
            tc   = sect_thick(r, c),
            b    = tpos_at(r),
            h    = H * pitch_fac(r),
            L    = sqrt(r * r + h * h),
            span = c / L * 180 / PI,     // angular chord span (deg)
            pc   = centre_off(r, c, L),
            z0   = rake_off(r))
        for (k = [0 : M - 1])
            let(j = k < K ? k : M - 1 - k,
                u = -cos(j * 180 / (K - 1)),          // cosine spacing, -1 = TE, +1 = LE
                n = k < K ? 0                          // flat pressure face
                          : edge_thickness + (tc - edge_thickness) * back_frac((1 - u) / 2, b))
            blade_pt(r, pc + u * span / 2, n, h, L, z0)
];

function blade_faces() = concat(
    // side surfaces (triangulated, since the helical quads are non-planar)
    [ for (i = [0 : N - 1]) for (j = [0 : M - 1])
        let(jn = (j + 1) % M,
            a = i * M + j, b = i * M + jn,
            c = (i + 1) * M + jn, d = (i + 1) * M + j)
        each [ [a, d, c], [a, c, b] ] ],
    // root cap
    [ for (j = [0 : K - 2]) each [
        [ j, j + 1, M - 1 - j ],
        [ j + 1, M - 2 - j, M - 1 - j ] ] ],
    // tip cap
    [ for (j = [0 : K - 2]) let(o = N * M) each [
        [ o + M - 1 - j, o + j + 1, o + j ],
        [ o + M - 1 - j, o + M - 2 - j, o + j + 1 ] ] ]
);

module blade() {
    polyhedron(points = blade_points(), faces = blade_faces(), convexity = 6);
}

// ---------------------------------------------------------------
// Hub and bore
// ---------------------------------------------------------------
module hub() {
    ch = hub_chamfer;
    translate([0, 0, -hub_len / 2])
        rotate_extrude()
            polygon([
                [0, 0],
                [hub_r_aft - ch, 0],
                [hub_r_aft, ch],
                [hub_r, hub_len - ch],
                [hub_r - ch, hub_len],
                [0, hub_len]
            ]);
}

module bore() {
    e = 1;
    len = hub_len + 2 * e;
    d_fwd = bore_d;
    d_aft = bore_taper > 0 ? bore_d - len / bore_taper : bore_d;
    translate([0, 0, -hub_len / 2 - e]) {
        cylinder(h = len, d1 = d_aft, d2 = d_fwd);
        if (keyway_width > 0 && keyway_depth > 0)
            translate([0, -keyway_width / 2, 0])
                cube([bore_d / 2 + keyway_depth, keyway_width, len]);
        if (shearpin_diameter > 0 && shearpin_length > 0)
            translate([0,0,len - 1 - shearpin_diameter/2]) {
                rotate([0,90,0])
                    cylinder(r=shearpin_diameter/2, h=shearpin_length, center=true);
                translate([-shearpin_length/2,-shearpin_diameter/2,0])
                    cube([shearpin_length, shearpin_diameter, shearpin_diameter]);
            }
    }
}

// ---------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------
module propeller() {
    difference() {
        union() {
            hub();
            for (b = [0 : num_blades - 1])
                rotate([0, 0, b * 360 / num_blades])
                    blade();
        }
        bore();
    }
}

if (hand == "left")
    mirror([1, 0, 0]) propeller();
else
    propeller();
