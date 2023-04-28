/*
 * Copyright (c) 2023-, Pedro Maciel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the FreeBSD Project.
 */


#include "CGAL/Delaunay_triangulation_on_sphere_2.h"
#include "CGAL/Delaunay_triangulation_on_sphere_traits_2.h"
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"


using K       = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits  = CGAL::Delaunay_triangulation_on_sphere_traits_2<K>;
using DToS2   = CGAL::Delaunay_triangulation_on_sphere_2<Traits>;
using Point_3 = Traits::Point_3;


int main(int argc, char* argv[]) {
    std::vector<Point_3> points{
        {2, 1, 1},   //
        {-2, 1, 1},  // not on the sphere
        {0, 1, 1},   //
        {1, 2, 1},   //
        {0, 1, 1},   // duplicate of #3
        {1, 0, 1},   //
        {1, 1, 2},   //
    };

    Traits traits({1, 1, 1}, 1);  // sphere center on (1,1,1), with radius 1

    DToS2 dtos(traits);
    for (const Point_3& pt : points) {
        std::cout << "Inserting (" << pt << ") at squared distance " << CGAL::squared_distance(pt, traits.center())
                  << " from the center of the sphere; is it on there sphere? "
                  << (traits.is_on_sphere(pt) ? "yes" : "no") << std::endl;
        dtos.insert(pt);

        std::cout << "After insertion, the dimension of the triangulation is: " << dtos.dimension() << "\n";
        std::cout << "It has:\n";
        std::cout << dtos.number_of_vertices() << " vertices\n";
        std::cout << dtos.number_of_edges() << " edges\n";
        std::cout << dtos.number_of_faces() << " solid faces\n";
        std::cout << dtos.number_of_ghost_faces() << " ghost faces\n" << std::endl;
    }

    CGAL::IO::write_OFF(std::cout, dtos, CGAL::parameters::stream_precision(17));

    return EXIT_SUCCESS;
}
