//
// Created by yuyttenh on 02/07/2021.
//

/**
 * @file geometry3d.h
 *
 * @brief Arbitrary exact and non-exact geometrical tests.
 */

#ifndef SWIFTSIM_GEOMETRY_3D_H
#define SWIFTSIM_GEOMETRY_3D_H

#include <gmp.h>
#include <math.h>

/**
 * @brief Auxiliary variables used by the arbirary exact tests. Since allocating
 * and deallocating these variables poses a significant overhead, they are best
 * reused.
 */
struct geometry3d {
  /*! @brief Arbitrary exact vertex coordinates */
  mpz_t aix, aiy, aiz, bix, biy, biz, cix, ciy, ciz, dix, diy, diz, eix, eiy,
      eiz;

  /*! @brief Temporary variables used to store relative vertex coordinates. */
  mpz_t s1x, s1y, s1z, s2x, s2y, s2z, s3x, s3y, s3z, s4x, s4y, s4z;

  /*! @brief Temporary variables used to store intermediate results. */
  mpz_t tmp1, tmp2, ab, bc, cd, da, ac, bd;

  /*! @brief Temporary variable used to store final exact results, before their
   *  sign is evaluated and returned as a finite precision integer. */
  mpz_t result;
};

/**
 * @brief Initialize the geometry3d object.
 *
 * This allocates and initialises the auxiliary arbitrary precision variables.
 *
 * @param g Geometry object.
 */
inline static void geometry3d_init(struct geometry3d* restrict g) {
  mpz_inits(g->aix, g->aiy, g->aiz, g->bix, g->biy, g->biz, g->cix, g->ciy,
            g->ciz, g->dix, g->diy, g->diz, g->eix, g->eiy, g->eiz, g->s1x,
            g->s1y, g->s1z, g->s2x, g->s2y, g->s2z, g->s3x, g->s3y, g->s3z,
            g->s4x, g->s4y, g->s4z, g->tmp1, g->tmp2, g->ab, g->bc, g->cd,
            g->da, g->ac, g->bd, g->result, NULL);
}

/**
 * @brief Deallocate all memory occupied by the geometry3d object.
 *
 * @param g Geometry object.
 */
inline static void geometry3d_destroy(struct geometry3d* restrict g) {
  mpz_clears(g->aix, g->aiy, g->aiz, g->bix, g->biy, g->biz, g->cix, g->ciy,
             g->ciz, g->dix, g->diy, g->diz, g->eix, g->eiy, g->eiz, g->s1x,
             g->s1y, g->s1z, g->s2x, g->s2y, g->s2z, g->s3x, g->s3y, g->s3z,
             g->s4x, g->s4y, g->s4z, g->tmp1, g->tmp2, g->ab, g->bc, g->cd,
             g->da, g->ac, g->bd, g->result, NULL);
}

inline static double geometry3d_orient(void) {
  // TODO
  return -1.;
}

/**
 * @brief Test the orientation of the tetrahedron that has the four given
 * points as vertex_indices.
 *
 * The test returns a positive result if the fourth vertex is below the plane
 * through the three other vertex_indices, with above the direction from which
 * the three points are ordered counterclockwise.
 *
 * E.g. if the four points are (0, 0, 0), (0, 0, 1), (0, 1, 0), and (1, 0, 0),
 * then this function returns 1.
 *
 * If the four points are exactly coplanar, then this function returns 0.
 *
 * @param a First vertex.
 * @param b Second vertex.
 * @param c Third vertex.
 * @param d Fourth vertex.
 * @return -1, 0, or 1, depending on the orientation of the tetrahedron.
 */
inline static int geometry3d_orient_exact(
    struct geometry3d* g, const unsigned long ax, const unsigned long ay,
    const unsigned long az, const unsigned long bx, const unsigned long by,
    const unsigned long bz, const unsigned long cx, const unsigned long cy,
    const unsigned long cz, const unsigned long dx, const unsigned long dy,
    const unsigned long dz) {

  /* store the input coordinates into the temporary large integer variables */
  mpz_set_ui(g->aix, ax);
  mpz_set_ui(g->aiy, ay);
  mpz_set_ui(g->aiz, az);

  mpz_set_ui(g->bix, bx);
  mpz_set_ui(g->biy, by);
  mpz_set_ui(g->biz, bz);

  mpz_set_ui(g->cix, cx);
  mpz_set_ui(g->ciy, cy);
  mpz_set_ui(g->ciz, cz);

  mpz_set_ui(g->dix, dx);
  mpz_set_ui(g->diy, dy);
  mpz_set_ui(g->diz, dz);

  /* compute large integer relative coordinates */
  mpz_sub(g->s1x, g->aix, g->dix);
  mpz_sub(g->s1y, g->aiy, g->diy);
  mpz_sub(g->s1z, g->aiz, g->diz);

  mpz_sub(g->s2x, g->bix, g->dix);
  mpz_sub(g->s2y, g->biy, g->diy);
  mpz_sub(g->s2z, g->biz, g->diz);

  mpz_sub(g->s3x, g->cix, g->dix);
  mpz_sub(g->s3y, g->ciy, g->diy);
  mpz_sub(g->s3z, g->ciz, g->diz);

  /* Compute the result in 3 steps */
  mpz_set_ui(g->result, 0);

  mpz_mul(g->tmp1, g->s2x, g->s3y);
  mpz_submul(g->tmp1, g->s3x, g->s2y);
  mpz_addmul(g->result, g->s1z, g->tmp1);

  mpz_mul(g->tmp1, g->s3x, g->s1y);
  mpz_submul(g->tmp1, g->s1x, g->s3y);
  mpz_addmul(g->result, g->s2z, g->tmp1);

  mpz_mul(g->tmp1, g->s1x, g->s2y);
  mpz_submul(g->tmp1, g->s2x, g->s1y);
  mpz_addmul(g->result, g->s3z, g->tmp1);

  return mpz_sgn(g->result);
}

inline static double geometry3d_in_sphere(void) {
  // TODO
  return -1.;
}

/**
 * @brief Check if the fifth given point is inside (-1) the circumsphere of the
 * tetrahedron formed by the other four given points.
 *
 * It is assumed that the first four points are the vertex_indices of a
 * positively oriented tetrahedron, as defined by a negative return value of
 * orient3d().
 *
 * If the fifth point is exactly on the circumsphere of the tetrahedron, this
 * functions returns 0.
 *
 * @param g Geometry struct
 * @param ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez Coordinates
 * of the vertex_indices of the tetrahedron and the test point (e).
 * @return -1, 0, or 1, depending on the outcome of the geometric test
 */
inline static int geometry3d_in_sphere_exact(
    struct geometry3d* restrict g, const unsigned long ax,
    const unsigned long ay, const unsigned long az, const unsigned long bx,
    const unsigned long by, const unsigned long bz, const unsigned long cx,
    const unsigned long cy, const unsigned long cz, const unsigned long dx,
    const unsigned long dy, const unsigned long dz, const unsigned long ex,
    const unsigned long ey, const unsigned long ez) {
  /* store the input coordinates into the temporary large integer variables */
  mpz_set_ui(g->aix, ax);
  mpz_set_ui(g->aiy, ay);
  mpz_set_ui(g->aiz, az);

  mpz_set_ui(g->bix, bx);
  mpz_set_ui(g->biy, by);
  mpz_set_ui(g->biz, bz);

  mpz_set_ui(g->cix, cx);
  mpz_set_ui(g->ciy, cy);
  mpz_set_ui(g->ciz, cz);

  mpz_set_ui(g->dix, dx);
  mpz_set_ui(g->diy, dy);
  mpz_set_ui(g->diz, dz);

  mpz_set_ui(g->eix, ex);
  mpz_set_ui(g->eiy, ey);
  mpz_set_ui(g->eiz, ez);

  /* compute large integer relative coordinates */
  mpz_sub(g->s1x, g->aix, g->eix);
  mpz_sub(g->s1y, g->aiy, g->eiy);
  mpz_sub(g->s1z, g->aiz, g->eiz);

  mpz_sub(g->s2x, g->bix, g->eix);
  mpz_sub(g->s2y, g->biy, g->eiy);
  mpz_sub(g->s2z, g->biz, g->eiz);

  mpz_sub(g->s3x, g->cix, g->eix);
  mpz_sub(g->s3y, g->ciy, g->eiy);
  mpz_sub(g->s3z, g->ciz, g->eiz);

  mpz_sub(g->s4x, g->dix, g->eix);
  mpz_sub(g->s4y, g->diy, g->eiy);
  mpz_sub(g->s4z, g->diz, g->eiz);

  /* compute intermediate values */
  mpz_mul(g->ab, g->s1x, g->s2y);
  mpz_submul(g->ab, g->s2x, g->s1y);

  mpz_mul(g->bc, g->s2x, g->s3y);
  mpz_submul(g->bc, g->s3x, g->s2y);

  mpz_mul(g->cd, g->s3x, g->s4y);
  mpz_submul(g->cd, g->s4x, g->s3y);

  mpz_mul(g->da, g->s4x, g->s1y);
  mpz_submul(g->da, g->s1x, g->s4y);

  mpz_mul(g->ac, g->s1x, g->s3y);
  mpz_submul(g->ac, g->s3x, g->s1y);

  mpz_mul(g->bd, g->s2x, g->s4y);
  mpz_submul(g->bd, g->s4x, g->s2y);

  /* compute the result in 4 steps */
  mpz_set_ui(g->result, 0);

  mpz_mul(g->tmp1, g->s4x, g->s4x);
  mpz_addmul(g->tmp1, g->s4y, g->s4y);
  mpz_addmul(g->tmp1, g->s4z, g->s4z);
  mpz_mul(g->tmp2, g->s1z, g->bc);
  mpz_submul(g->tmp2, g->s2z, g->ac);
  mpz_addmul(g->tmp2, g->s3z, g->ab);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s3x, g->s3x);
  mpz_addmul(g->tmp1, g->s3y, g->s3y);
  mpz_addmul(g->tmp1, g->s3z, g->s3z);
  mpz_mul(g->tmp2, g->s4z, g->ab);
  mpz_addmul(g->tmp2, g->s1z, g->bd);
  mpz_addmul(g->tmp2, g->s2z, g->da);
  mpz_submul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s2x, g->s2x);
  mpz_addmul(g->tmp1, g->s2y, g->s2y);
  mpz_addmul(g->tmp1, g->s2z, g->s2z);
  mpz_mul(g->tmp2, g->s3z, g->da);
  mpz_addmul(g->tmp2, g->s4z, g->ac);
  mpz_addmul(g->tmp2, g->s1z, g->cd);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s1x, g->s1x);
  mpz_addmul(g->tmp1, g->s1y, g->s1y);
  mpz_addmul(g->tmp1, g->s1z, g->s1z);
  mpz_mul(g->tmp2, g->s2z, g->cd);
  mpz_submul(g->tmp2, g->s3z, g->bd);
  mpz_addmul(g->tmp2, g->s4z, g->bc);
  mpz_submul(g->result, g->tmp1, g->tmp2);

  return mpz_sgn(g->result);
}

/**
 * @brief Compute the coordinates of the circumcenter of the tetrahedron
 * (v0, v1, v2, v3).
 *
 * See https://mathworld.wolfram.com/Circumsphere.html
 *
 * @param v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z Coordinates
 * of the corners of the tetrahedron.
 * @param circumcenter (Returned) coordinates of center of circumsphere
 * @param Ry (Returned) y coordinate of center of circumsphere
 * @param Rz (Returned) z coordinate of center of circumsphere
 */
static inline void geometry3d_compute_circumcenter(
    double v0x, double v0y, double v0z, double v1x, double v1y, double v1z,
    double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
    double* circumcenter) {
  /* Compute relative coordinates */
  const double r1x = v1x - v0x;
  const double r1y = v1y - v0y;
  const double r1z = v1z - v0z;
  const double r2x = v2x - v0x;
  const double r2y = v2y - v0y;
  const double r2z = v2z - v0z;
  const double r3x = v3x - v0x;
  const double r3y = v3y - v0y;
  const double r3z = v3z - v0z;

  /* Compute squared norm of relative coordinates */
  const double r1_sqrd = r1x * r1x + r1y * r1y + r1z * r1z;
  const double r2_sqrd = r2x * r2x + r2y * r2y + r2z * r2z;
  const double r3_sqrd = r3x * r3x + r3y * r3y + r3z * r3z;

  const double Dx = r1_sqrd * (r2y * r3z - r3y * r2z) -
                    r2_sqrd * (r1y * r3z - r3y * r1z) +
                    r3_sqrd * (r1y * r2z - r2y * r1z);
  const double Dy = -r1_sqrd * (r2x * r3z - r3x * r2z) +
                    r2_sqrd * (r1x * r3z - r3x * r1z) -
                    r3_sqrd * (r1x * r2z - r2x * r1z);
  const double Dz = r1_sqrd * (r2x * r3y - r3x * r2y) -
                    r2_sqrd * (r1x * r3y - r3x * r1y) +
                    r3_sqrd * (r1x * r2y - r2x * r1y);

  const double a = r1x * (r2y * r3z - r3y * r2z) -
                   r2x * (r1y * r3z - r3y * r1z) +
                   r3x * (r1y * r2z - r2y * r1z);

  const double denominator = 2. * a;
  circumcenter[0] = Dx / denominator + v0x;
  circumcenter[1] = Dy / denominator + v0y;
  circumcenter[2] = Dz / denominator + v0z;
}

inline static double geometry3d_compute_area_triangle(double ax, double ay,
                                                      double az, double bx,
                                                      double by, double bz,
                                                      double cx, double cy,
                                                      double cz) {
  const double abx = bx - ax;
  const double aby = by - ay;
  const double abz = bz - az;
  const double acx = cx - ax;
  const double acy = cy - ay;
  const double acz = cz - az;

  const double Dx = aby * acz - abz * acy;
  const double Dy = abz * acx - abx * acz;
  const double Dz = abx * acy - aby * acx;

  return sqrt(Dx * Dx + Dy * Dy + Dz * Dz) / 2.;
}

inline static void geometry3d_compute_centroid_triangle(
    double ax, double ay, double az, double bx, double by, double bz, double cx,
    double cy, double cz, double* result) {
  result[0] = (ax + bx + cx) / 3.;
  result[1] = (ay + by + cy) / 3.;
  result[2] = (az + bz + cz) / 3.;
}

inline static double geometry3d_compute_volume_tetrahedron(
    double ax, double ay, double az, double bx, double by, double bz, double cx,
    double cy, double cz, double dx, double dy, double dz) {
  /* Compute relative coordinates */
  const double dax = ax - dx;
  const double day = ay - dy;
  const double daz = az - dz;
  const double dbx = bx - dx;
  const double dby = by - dy;
  const double dbz = bz - dz;
  const double dcx = cx - dx;
  const double dcy = cy - dy;
  const double dcz = cz - dz;

  /* compute (b - d) x (c - d) */
  const double cross_x = dby * dcz - dcy * dbz;
  const double cross_y = dbx * dcz - dcx * dbz;
  const double cross_z = dbx * dcy - dcx * dby;

  return fabs(dax * cross_x - day * cross_y + daz * cross_z) / 6.;
}

inline static void geometry3d_compute_centroid_tetrahedron(
    double ax, double ay, double az, double bx, double by, double bz, double cx,
    double cy, double cz, double dx, double dy, double dz, double* result) {
  result[0] = (ax + bx + cx + dx) / 4.;
  result[1] = (ay + by + cy + dy) / 4.;
  result[2] = (az + bz + cz + dz) / 4.;
}

inline static double geometry3d_compute_centroid_volume_tetrahedron(
    double ax, double ay, double az, double bx, double by, double bz, double cx,
    double cy, double cz, double dx, double dy, double dz, double* result) {
  geometry3d_compute_centroid_tetrahedron(ax, ay, az, bx, by, bz, cx, cy, cz,
                                          dx, dy, dz, result);
  return geometry3d_compute_volume_tetrahedron(ax, ay, az, bx, by, bz, cx, cy,
                                               cz, dx, dy, dz);
}

inline static double geometry3d_compute_centroid_area(
    const double* restrict points, int n_points, double* result) {
  result[0] = 0.;
  result[1] = 0.;
  result[2] = 0.;
  for (int i = 0; i < n_points; i++) {
    result[0] += points[3 * i];
    result[1] += points[3 * i + 1];
    result[2] += points[3 * i + 2];
  }

  if (n_points < 2) return 0.;

  double V = 0.;
  const double v0x = points[0];
  const double v0y = points[1];
  const double v0z = points[2];
  for (int i = 2; i < n_points; i++) {
    V += geometry3d_compute_area_triangle(
        v0x, v0y, v0z, points[3 * i - 3], points[3 * i - 2], points[3 * i - 1],
        points[3 * i], points[3 * i + 1], points[3 * i + 2]);
  }
  return V;
}


#endif  // SWIFTSIM_GEOMETRY_3D_H
