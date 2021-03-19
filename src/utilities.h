/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef SWIFT_UTILITIES_H
#define SWIFT_UTILITIES_H

/**
 * @brief Search for a value in a monotonically increasing array to find the
 *      index such that array[index] < value < array[index + 1]
 *
 * @param x The value to find
 * @param array The array to search
 * @param n The length of the array
 *
 * Return -1 and n for x below and above the array edge values respectively.
 */
INLINE static int find_value_in_monot_incr_array(const float x,
                                                 const float *array,
                                                 const int n) {

  int index_mid, index_low = 0, index_high = n;

  // Until array[index_low] < x < array[index_high=index_low+1]
  while (index_high - index_low > 1) {
    index_mid = (index_high + index_low) / 2;  // Middle index

    // Replace the low or high index with the middle
    if (array[index_mid] <= x)
      index_low = index_mid;
    else
      index_high = index_mid;
  }

  // Set index with the found index_low or an error value if outside the array
  if (x < array[0])
    return -1;
  else if (array[n - 1] <= x)
    return n;
  else
    return index_low;
}

INLINE static float imbalance_statistic_q(const float N) {

  // montecarlo N = 1M, N_neig 5 to 50
  /* quantile 99 */
  /*const float a = -1.53377904f;
  const float b = 0.23923416f;
  const float c = 1.5047929f;*/

  /* quantile 95 */ /*
  const float a = -0.31934966f;
  const float b = 0.31051137f;
  const float c = 1.24962448f;*/

  /* quantile 90 */
  /*const float a = 0.04169535f;
  const float b = 0.33453515f;
  const float c = 1.1176431f;*/

  /* quantile 80 */
  /*const float a = 0.3720982f;
  const float b = 0.32889313f;
  const float c = 0.96356348f;*/

  /* quantile 70 */
  const float a = 0.69223044f;
  const float b = 0.28108532f;
  const float c = 0.85737999f;


  const float N_inv = 1.f / N;
  float q99 = 0.f;

  q99 = a*N_inv*N_inv + b*N_inv + c;

  return q99;

}



#endif /* SWIFT_UTILITIES_H */
