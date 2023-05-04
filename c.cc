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


#include <array>
#include <fstream>
#include <random>
#include <vector>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>


struct Point : std::array<double, 3> {
    Point(double x, double y, double z) : array{x, y, z} {}

    double& x = operator[](0);
    double& y = operator[](1);
    double& z = operator[](2);
};


void c() {
    auto rnd = []() -> Point {
        static std::random_device rd{};
        static std::mt19937 gen{rd()};
        static std::normal_distribution<double> d{0., 2};

        for (;;) {
            auto x = d(gen);
            auto y = d(gen);
            auto z = d(gen);

            if (auto norm = std::sqrt(x * x + y * y + z * z); norm > 0) {
                return {x / norm, y / norm, z / norm};
            }
        }
    };

    size_t N = 100;

    std::vector<coordT> coords;
    coords.reserve(3 * N);

    for (size_t c = 0; c < N; c++) {
        auto p = rnd();
        coords.emplace_back(p.x);
        coords.emplace_back(p.y);
        coords.emplace_back(p.z);
    }

    orgQhull::Qhull qh;
    qh.runQhull("", 3, static_cast<int>(N), coords.data(), "QJ");

#if 0
    auto mm =
        std::minmax(qh.beginVertex(), qh.endVertex(),
                    [](const orgQhull::QhullVertex& a, const orgQhull::QhullVertex& b) { return a.id() < b.id(); });
    std::cout << "min: " << mm.first.id() << " max: " << mm.second.id() << std::endl;
#endif

    {
        std::ofstream out("trish.msh");

        out << "$MeshFormat\n"
            << "2.2 0 " << sizeof(double) << '\n'
            << "$EndMeshFormat\n";

        out << "$Nodes\n" << qh.vertexCount() << '\n';
        size_t i = 0;
        for (const auto& v : qh.vertexList()) {
            out << v.id() << ' ' << v.point();
        }
        out << "$EndNodes\n";

        out << "$Elements\n" << qh.facetCount() << '\n';
        size_t j = 0;
        for (const auto& facet : qh.facetList()) {
            assert(facet.isSimplicial());
            out << j++ << " 2 4 1 1 1 0 " << (facet.vertices().at(0).id()) << ' ' << (facet.vertices().at(1).id())
                << ' ' << (facet.vertices().at(2).id()) << '\n';
        }
        out << "$EndElements\n";
    }
}
