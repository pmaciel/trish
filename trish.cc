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


#include <unordered_map>

#include "CGAL/Delaunay_triangulation_on_sphere_2.h"
#include "CGAL/Delaunay_triangulation_on_sphere_traits_2.h"
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"


using Kernel        = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits        = CGAL::Delaunay_triangulation_on_sphere_traits_2<Kernel>;
using Triangulation = CGAL::Delaunay_triangulation_on_sphere_2<Traits>;


int main(int argc, char* argv[]) {
    std::vector<Traits::Point_3> points{
        {2, 1, 1},   //
        {-2, 1, 1},  // not on the sphere
        {0, 1, 1},   //
        {1, 2, 1},   //
        {0, 1, 1},   // duplicate of #3
        {1, 0, 1},   //
        {1, 1, 2},   //
    };

    Traits traits({1, 1, 1}, 1);  // sphere center on (1,1,1), with radius 1

    Triangulation tri(points.begin(), points.end(), traits);

    std::cout << "dimension: " << tri.dimension() << '\n'
              << "number_of_vertices: " << tri.number_of_vertices() << '\n'
              << "number_of_edges: " << tri.number_of_edges() << '\n'
              << "number_of_faces: " << tri.number_of_faces() << '\n'
              << "number_of_ghost_faces: " << tri.number_of_ghost_faces() << std::endl;


    std::unordered_map<Triangulation::Vertex_handle, Triangulation::size_type> index;
    Triangulation::size_type i = 0;
    for (auto it = tri.vertices_begin(); it != tri.vertices_end(); ++it) {
        index[it] = i++;
    }

    CGAL_triangulation_assertion(i == tri.number_of_vertices());

    for (auto it = tri.all_faces_begin(); it != tri.all_faces_end(); ++it) {
        if (!it->is_ghost()) {
            std::cout << index[it->vertex(0)] << " " << index[it->vertex(1)] << " " << index[it->vertex(2)] << '\n';
        }
    }


    return 0;
}
