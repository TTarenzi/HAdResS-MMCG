/*
<<<<<<< HEAD
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
=======
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 */

#ifndef _GMX_RANDOM_H_
#define _GMX_RANDOM_H_
<<<<<<< HEAD

=======
#include "visibility.h"
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
#include <stdio.h>
#include "types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Abstract datatype for a random number generator
 *
<<<<<<< HEAD
 * This is a handle to the full state of a random number generator. 
=======
 * This is a handle to the full state of a random number generator.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 * You can not access anything inside the gmx_rng structure outside this
 * file.
 */
typedef struct gmx_rng *
<<<<<<< HEAD
gmx_rng_t;
=======
    gmx_rng_t;
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/*! \brief Returns the size of the RNG integer data structure
 *
 * Returns the size of the RNG integer data structure.
 * \threadsafe Yes.
 */
<<<<<<< HEAD
=======
GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
int
gmx_rng_n(void);


/*! \brief Create a new RNG, seeded from a single integer.
 *
 * If you dont want to pick a seed, just call it as
 * rng=gmx_rng_init(gmx_rng_make_seed()) to seed it from
 * the system time or a random device.
 *
 * \param seed Random seed, unsigned 32-bit integer.
 *
<<<<<<< HEAD
 * \return Reference to a random number generator, or NULL if there was an 
=======
 * \return Reference to a random number generator, or NULL if there was an
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *         error.
 *
 * \threadsafe Yes.
 */
<<<<<<< HEAD
gmx_rng_t 
=======
GMX_LIBGMX_EXPORT
gmx_rng_t
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
gmx_rng_init(unsigned int seed);


/*! \brief Generate a 'random' RNG seed.
 *
 * This routine tries to get a seed from /dev/random if present,
 * and if not it uses time-of-day and process id to generate one.
 *
<<<<<<< HEAD
 * \return 32-bit unsigned integer random seed. 
=======
 * \return 32-bit unsigned integer random seed.
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *
 * Tip: If you use this in your code, it is a good idea to write the
 * returned random seed to a logfile, so you can recreate the exact sequence
 * of random number if you need to reproduce your run later for one reason
 * or another.
 *
 * \threadsafe Yes.
 */
<<<<<<< HEAD
=======
GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
unsigned int
gmx_rng_make_seed(void);


/*! \brief Initialize a RNG with 624 integers (>32 bits of entropy).
 *
 *  The Mersenne twister RNG used in Gromacs has an extremely long period,
 *  but when you only initialize it with a 32-bit integer there are only
 *  2^32 different possible sequences of number - much less than the generator
 *  is capable of.
 *
 *  If you really need the full entropy, this routine makes it possible to
<<<<<<< HEAD
 *  initialize the RNG with up to 624 32-bit integers, which will give you 
=======
 *  initialize the RNG with up to 624 32-bit integers, which will give you
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *  up to 2^19968 bits of entropy.
 *
 *  \param seed Array of unsigned integers to form a seed
 *  \param seed_length Number of integers in the array, up to 624 are used.
 *
<<<<<<< HEAD
 * \return Reference to a random number generator, or NULL if there was an 
=======
 * \return Reference to a random number generator, or NULL if there was an
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *         error.
 *
 * \threadsafe Yes.
 */
<<<<<<< HEAD
gmx_rng_t 
gmx_rng_init_array(unsigned int    seed[], 
		   int             seed_length);
=======
gmx_rng_t
gmx_rng_init_array(unsigned int    seed[],
                   int             seed_length);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/*! \brief Release resources of a RNG
 *
 *  This routine destroys a random number generator and releases all
 *  resources allocated by it.
 *
 *  \param rng Handle to random number generator previously returned by
 *		       gmx_rng_init() or gmx_rng_init_array().
 *
<<<<<<< HEAD
 * \threadsafe Function itself is threadsafe, but you should only destroy a 
 *             certain RNG once (i.e. from one thread).
 */
=======
 * \threadsafe Function itself is threadsafe, but you should only destroy a
 *             certain RNG once (i.e. from one thread).
 */
GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
void
gmx_rng_destroy(gmx_rng_t rng);


/*! \brief Get the state of a RNG
 *
 *  This routine stores the random state in mt and mti, mt should have
 *  a size of at least 624, mt of 1.
 *
 *  \param rng Handle to random number generator previously returned by
 *		       gmx_rng_init() or gmx_rng_init_array().
 */
<<<<<<< HEAD
void
gmx_rng_get_state(gmx_rng_t rng, unsigned int *mt,int *mti);
=======
GMX_LIBGMX_EXPORT
void
gmx_rng_get_state(gmx_rng_t rng, unsigned int *mt, int *mti);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/*! \brief Set the state of a RNG
 *
 *  This routine sets the random state from mt and mti, mt should have
 *  a size of at least 624.
 *
 *  \param rng Handle to random number generator previously returned by
 *		       gmx_rng_init() or gmx_rng_init_array().
 */
<<<<<<< HEAD
void
gmx_rng_set_state(gmx_rng_t rng, unsigned int *mt,int mti);
=======
GMX_LIBGMX_EXPORT
void
gmx_rng_set_state(gmx_rng_t rng, unsigned int *mt, int mti);
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2


/*! \brief Random 32-bit integer from a uniform distribution
 *
 *  This routine returns a random integer from the random number generator
 *  provided, and updates the state of that RNG.
 *
 *  \param rng Handle to random number generator previously returned by
 *		       gmx_rng_init() or gmx_rng_init_array().
 *
 *  \return 32-bit unsigned integer from a uniform distribution.
 *
<<<<<<< HEAD
 *  \threadsafe Function yes, input data no. You should not call this function 
 *	        from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle 
 *              with a mutex every time we need a random number - that would 
 *              slow the routine down a factor 2-5. There are two simple 
 *		solutions: either use a mutex and lock it before calling
 *              the function, or use a separate RNG handle for each thread.
 */
unsigned int 
=======
 *  \threadsafe Function yes, input data no. You should not call this function
 *	        from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle
 *              with a mutex every time we need a random number - that would
 *              slow the routine down a factor 2-5. There are two simple
 *		solutions: either use a mutex and lock it before calling
 *              the function, or use a separate RNG handle for each thread.
 */
GMX_LIBGMX_EXPORT
unsigned int
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
gmx_rng_uniform_uint32(gmx_rng_t rng);


/*! \brief Random gmx_real_t 0<=x<1 from a uniform distribution
 *
<<<<<<< HEAD
 *  This routine returns a random floating-point number from the 
=======
 *  This routine returns a random floating-point number from the
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *  random number generator provided, and updates the state of that RNG.
 *
 *  \param rng Handle to random number generator previously returned by
 *		       gmx_rng_init() or gmx_rng_init_array().
 *
 *  \return floating-point number 0<=x<1 from a uniform distribution.
 *
<<<<<<< HEAD
 *  \threadsafe Function yes, input data no. You should not call this function 
 *		from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle 
 *              with a mutex every time we need a random number - that would 
 *              slow the routine down a factor 2-5. There are two simple 
 *		solutions: either use a mutex and lock it before calling
 *              the function, or use a separate RNG handle for each thread.
 */
=======
 *  \threadsafe Function yes, input data no. You should not call this function
 *		from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle
 *              with a mutex every time we need a random number - that would
 *              slow the routine down a factor 2-5. There are two simple
 *		solutions: either use a mutex and lock it before calling
 *              the function, or use a separate RNG handle for each thread.
 */
GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
real
gmx_rng_uniform_real(gmx_rng_t rng);


/*! \brief Random gmx_real_t from a gaussian distribution
 *
<<<<<<< HEAD
 *  This routine returns a random floating-point number from the 
 *  random number generator provided, and updates the state of that RNG.
 *  
=======
 *  This routine returns a random floating-point number from the
 *  random number generator provided, and updates the state of that RNG.
 *
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *  The Box-Muller algorithm is used to provide gaussian random numbers. This
 *  is not the fastest known algorithm for gaussian numbers, but in contrast
 *  to the alternatives it is very well studied and you can trust the returned
 *  random numbers to have good properties and no correlations.
 *
 *  \param rng Handle to random number generator previously returned by
 *			  gmx_rng_init() or gmx_rng_init_array().
 *
<<<<<<< HEAD
 *  \return Gaussian random floating-point number with average 0.0 and 
=======
 *  \return Gaussian random floating-point number with average 0.0 and
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *	    standard deviation 1.0. You can get any average/mean you want
 *          by first multiplying with the desired average and then adding
 *          the average you want.
 *
<<<<<<< HEAD
 *  \threadsafe Function yes, input data no. You should not call this function 
 *		from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle 
 *              with a mutex every time we need a random number - that would 
 *              slow the routine down a factor 2-5. There are two simple 
=======
 *  \threadsafe Function yes, input data no. You should not call this function
 *		from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle
 *              with a mutex every time we need a random number - that would
 *              slow the routine down a factor 2-5. There are two simple
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *		solutions: either use a mutex and lock it before calling
 *              the function, or use a separate RNG handle for each thread.
 *
 *  It works perfectly to mix calls to get uniform and gaussian random numbers
 *  from the same generator, but since it will affect the sequence of returned
 *  numbers it is probably better to use separate random number generator
 *  structures.
 */
<<<<<<< HEAD
=======
GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
real
gmx_rng_gaussian_real(gmx_rng_t rng);



/* Return a new gaussian random number with expectation value
 * 0.0 and standard deviation 1.0. This routine uses a table
 * lookup for maximum speed.
 *
 * WARNING: The lookup table is 16k by default, which means
 *          the granularity of the random numbers is coarser
 *	    than what you get from gmx_rng_gauss_real().
 *          In most cases this is no problem whatsoever,
 *          and it is particularly true for BD/SD integration.
<<<<<<< HEAD
 *	    Note that you will NEVER get any really extreme 
=======
 *	    Note that you will NEVER get any really extreme
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
 *          numbers: the maximum absolute value returned is
 *          4.0255485.
 *
 * threadsafe: yes
 */
<<<<<<< HEAD
=======
GMX_LIBGMX_EXPORT
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
real
gmx_rng_gaussian_table(gmx_rng_t rng);

#ifdef __cplusplus
}
#endif

#endif /* _GMX_RANDOM_H_ */
<<<<<<< HEAD

=======
>>>>>>> 2aa1ce7f6f65a45829ff249723a2667c864a68b2
