/******************************************************************************
 *  PlaTec, a 2D terrain generator based on plate tectonics
 *  Copyright (C) 2012- Lauri Viitanen
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, see http://www.gnu.org/licenses/
 *****************************************************************************/

/** @file	hmapgen_sqrdmd.c
 *  @brief	Contains functions to generate fractal height maps.
 *
 *  @author Lauri Viitanen
 *  @date 2011-08-09
 */
#include <stdlib.h>
#include <cstdio>

#include "sqrdmd.hpp"

#define CALC_SUM(a, b, c, d)\
{\
	sum = ((a) + (b) + (c) + (d)) * 0.25f;\
	sum = sum + slope * ((rand() << 1) - RAND_MAX);\
}

#define SAVE_SUM(a)\
{\
	bool isZero = (int)map[a] == 0;  \
	if (isZero) { \
		map[a] = sum; \
	} \
}

void normalize(float* arr, int size)
{
	float min = arr[0], max = arr[0], diff;

	for (int i = 1; i < size; ++i)
	{
		min = min < arr[i] ? min : arr[i];
		max = max > arr[i] ? max : arr[i];
	}

	diff = max - min;

	if (diff > 0)
		for (int i = 0; i < size; ++i)
			arr[i] = (arr[i] - min) / diff;
}

int sqrdmd(float* map, const int size, float rgh)
{
	return sqrdmd(map, size, size, rgh);
}

int sqrdmd(float* map, const int width, const int height, float rgh)
{	
	const int size = width; // temporary: it should be removed later
	const int full_size = width * height;

	int i;
	int x, y, dx, dy;
	int x0, x1, y0, y1;
	int p0, p1, p2, p3;
	int step; // temporary: it should be removed later
	int step_x, step_y, line_jump, masked;
	float slope, sum, center_sum;

	i = 0;
	slope  = rgh;
	step   = size   - 1;
	step_x = width  - 1;
	step_y = height - 1;

	/* Calculate midpoint ("diamond step"). */
	dy = step_y * width; // start of last row
	CALC_SUM(map[0], map[step_x], map[dy], map[dy + step_x]);
	SAVE_SUM(i);

	center_sum = sum;

	/* Calculate each sub diamonds' center points ("square step"). */

	/* Top row. */
	p0 = step_x >> 1;
	CALC_SUM(map[0], map[step_x], center_sum, center_sum);
	SAVE_SUM(p0);

	/* Left column. */
	p1 = (step_y >> 1) * width;
	CALC_SUM(map[0], map[dy], center_sum, center_sum);
	SAVE_SUM(p1);

	map[full_size + p0 - width] = map[p0]; /* Copy top val into btm row. */
	map[p1 + width - 1] = map[p1]; /* Copy left value into right column. */

	slope *= rgh;
	step   >>= 1; // temporary
	step_x >>= 1;
	step_y >>= 1;

	while (step_x > 1 && step_y > 1)  /* Enter the main loop. */
	{
		/*************************************************************
		 * Calc midpoint of sub squares on the map ("diamond step"). *
		 *************************************************************/

		dx = step_x;
		dy = step_y * width;
		
		i  = (step_y >> 1) * (width ) + (step_x >> 1);		
		line_jump = step_y * width + 1 + step_x - width;

		for (y0 = 0, y1 = dy; y1 < width * height; y0 += dy, y1 += dy)
		{
			for (x0 = 0, x1 = dx; x1 < width; x0 += dx, x1 += dx, i += step_x)
			{
                CALC_SUM(map[y0 + x0], map[y0 + x1], map[y1 + x0], map[y1 + x1]);
				SAVE_SUM(i);
			}

			/* There's additional step taken at the end of last
			 * valid loop. That step actually isn't valid because
			 * the row ends right then. Thus we are forced to
			 * manually remove it after the loop so that 'i'
			 * points again to the index accessed last.
			 */
			i += line_jump - step_x;
		}

		/**************************************************************
		 * Calculate each sub diamonds' center point ("square step").
		 * Diamond gets its left and right vertices from the square
		 * corners of last iteration and its top and bottom vertices
		 * from the "diamond step" we just performed.					 
		 *************************************************************/

		i  = step >> 1; // temporary: should be removed
		int ix = step_x / 2, iy = step_y / 2;
		p0 = step_x;  /* right */
		p1 = iy * width + ix;  /* bottom */
		p2 = 0;  /* left */
		p3 = full_size + ix - (iy + 1) * width; /* top (wrapping edges) */

		/* Calculate "diamond" values for top row in map. */
		while (p0 < width)
		{
            CALC_SUM(map[p0], map[p1], map[p2], map[p3]);
			SAVE_SUM(i);
			/* Copy it into bottom row. */
			map[full_size + ix - width] = map[ix];

			p0 += step_x; 
			p1 += step_x; 
			p2 += step_x;
			p3 += step_x; 
			i  += step_x;
			ix += step_x;
			iy += step_y;
		}

		/* Now that top row's values are calculated starting from
		 * 'y = step >> 1' both saves us from recalculating same things
		 * twice and guarantees that data will not be read beyond top
		 * row of map. 'size - (step >> 1)' guarantees that data will
		 * not be read beyond bottom row of map.
		 */
		int temp;
		for (y = step_y >> 1, temp = 0; 
			 y < height - (step_y >> 1);
			 y += step_y >> 1, temp = !temp)
		{
			int i;
			p0 = step_x >> 1;  /* right */
			p1 = (step_y >> 1) * width;  /* bottom */
			p2 = -p0;  /* left */
			p3 = -p1;  /* top */

			/* For even rows add step/2. Otherwise add nothing. */
			x = i = p0 * temp;  /* Init 'x' while it's easy. */
			i += y * width;  /* Move 'i' into correct row. */

			p0 += i;
			p1 += i;
			/* For odd rows p2 (left) wraps around map edges. */
			p2 += i + (width - 1) * !temp;  
			p3 += i;

			/* size - (step >> 1) guarantees that data will not be
			 * read beyond rightmost column of map. */
			for (; x < width - (step_x >> 1); x += step_x)
			{
				CALC_SUM(map[p0], map[p1], map[p2], map[p3]);
				SAVE_SUM(i);

				p0 += step;
				p1 += step;
				p2 += step;
				p3 += step;
				i  += step;

				/* if we start from leftmost column -> left
				 * point (p2) is going over the right border ->
				 * wrap it around into the beginning of
				 * previous rows left line. */
				p2 -= (size - 1) * !x;
			}

			/* copy rows first element into its last */
			i = y * size;
			map[i + size - 1] = map[i];
		}

		slope *= rgh; /* reduce amount of randomness for next round */
		step   /= 2; /* split squares and diamonds in half */
		step_x /= 2;
		step_y /= 2;
	}

	return (0);
}
