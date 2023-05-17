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


#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>


struct Point : std::array<double, 2> {
    Point(value_type x, value_type y) : array{x, y} {}
    Point(const Point& P) : array(P) {}
    Point(Point&& P) : array(P) {}

    Point& operator=(const Point& P) {
        x = P.x;
        y = P.y;
        return *this;
    }

    Point& operator=(Point&& P) {
        x = P.x;
        y = P.y;
        return *this;
    }

    value_type& x = operator[](0);
    value_type& y = operator[](1);

    friend Point operator+(const Point& a, const Point& b) { return {a.x + b.x, a.y + b.y}; }
    friend Point operator-(const Point& a, const Point& b) { return {a.x - b.x, a.y - b.y}; }

    friend std::ostream& operator<<(std::ostream& out, const Point& P) {
        return out << '{' << P.x << ", " << P.y << '}';
    }
};


struct Triangle;

using TrianglePtr = std::shared_ptr<Triangle>;


struct Edge {
    Edge(int v1, int v2, TrianglePtr _tri = nullptr) : v{v1, v2}, tri(_tri) {}
    const std::array<int, 3> v;  // Vertices
    TrianglePtr tri;

    bool operator==(const Edge& other) const { return v == other.v && tri == other.tri; }
};


std::vector<Point> __points;


struct Triangle {
    Triangle(int v1, int v2, int v3) : v{v1, v2, v3}, circumcenter_(__points[v1]) {
        n.fill(nullptr);

        auto AB = __points[v2] - __points[v1];
        auto AC = __points[v3] - __points[v1];
        auto d  = 2. * (AB.x * AC.y - AB.y * AC.x);

        auto b = AB.x * AB.x + AB.y * AB.y;
        auto c = AC.x * AC.x + AC.y * AC.y;

        auto x = (AC.y * b - AB.y * c) / d;
        auto y = (AB.x * c - AC.x * b) / d;

        circumcenter_.x += x;
        circumcenter_.y += y;

        circumradius2_ = x * x + y * y;
        circumradius_  = std::sqrt(circumradius2_);
    }

    const Point& A() const { return __points[v[0]]; }
    const Point& B() const { return __points[v[1]]; }
    const Point& C() const { return __points[v[2]]; }

    const Point& circumcenter() const { return circumcenter_; }
    Point::value_type circumradius() const { return circumradius_; }
    Point::value_type circumradius2() const { return circumradius2_; }

    bool in_circumcircle(const Point& P) {
        const auto D = P - circumcenter_;
        return D.x * D.x + D.y * D.y <= circumradius2_;
    }

    void set_edge(const Edge edge, const TrianglePtr T) {
        // Set the edge neighbour that matches "edge" to T

        for (int i : {0, 1, 2}) {
            if (edge.v[0] == v[i] && edge.v[1] == v[(i + 1) % 3]) {
                n[(i + 2) % 3] = T;
                return;
            }
        }
    }

    const std::array<int, 3> v;    // Vertices
    std::array<TrianglePtr, 3> n;  // Neighbours

private:
    Point circumcenter_;
    Point::value_type circumradius_;
    Point::value_type circumradius2_;

    friend std::ostream& operator<<(std::ostream& out, const Triangle& T) {
        return out << "{A: " << T.A() << ", B: " << T.B() << ", C: " << T.C() << ", circumcenter: " << T.circumcenter_
                   << ", circumradius: " << T.circumradius_ << '}';
    }
};


class Triangulation {
public:
    std::vector<TrianglePtr> triangles;

    explicit Triangulation(double R) {
        __points.emplace_back(-R, -R);
        __points.emplace_back(R, -R);
        __points.emplace_back(R, R);
        __points.emplace_back(-R, R);

        // Form the frame
        auto T1 = std::make_shared<Triangle>(0, 3, 1);
        auto T2 = std::make_shared<Triangle>(2, 1, 3);

        T1->n[0] = T2;
        T2->n[0] = T1;

        triangles.push_back(T1);
        triangles.push_back(T2);
    }

    void insert(Point&& P) {
        __points.push_back(P);
        auto p = static_cast<int>(__points.size() - 1);

        // Find bad triangles and their boundary
        std::vector<TrianglePtr> bad;
        std::vector<Edge> boundary;

        // FIXME: could be improved
        for (const auto& T : triangles) {
            if (T->in_circumcircle(P)) {
                bad.push_back(T);
            }
        }

        {
            auto T = bad[0];
            int e  = 0;

            while (true) {
                if (!boundary.empty() && boundary.front() == boundary.back()) {
                    break;
                }

                if (std::find(bad.begin(), bad.end(), T->n[e]) != bad.end()) {
                    // If this edge is shared with a triangle in bad_triangles, set current triangle
                    auto last = T;
                    T         = T->n[e];

                    auto p = static_cast<int>(std::distance(T->n.begin(), std::find(T->n.begin(), T->n.end(), last)));
                    e      = (p + 1) % 3;
                }
                else {
                    // Found an edge that is on the boundary, add to list
                    boundary.emplace_back(T->v[(e + 1) % 3], T->v[(e + 2) % 3], T->n[e]);
                    e = (e + 1) % 3;
                }
            }

            boundary.pop_back();
        }


        // Remove bad triangles
        for (const auto& T : bad) {
            triangles.erase(std::remove(triangles.begin(), triangles.end(), T), triangles.end());
        }

        // Retriangle the hole just created
        std::vector<TrianglePtr> more;
        for (const auto& edge : boundary) {
            auto a = edge.v[0];
            auto b = edge.v[1];

            auto T = std::make_shared<Triangle>(p, a, b);

            T->n[0] = edge.tri;  // To neighbour

            if (edge.tri) {
                T->n[0]->set_edge(Edge(b, a, nullptr), T);  // From neighbour
            }

            more.push_back(T);
        }

        // Link the new triangles to each other
        auto N = static_cast<int>(more.size());
        for (int i = 0; i < N; i++) {
            more[i]->n[2] = more[((i - 1) % N + N) % N];
            more[i]->n[1] = more[(i + 1) % N];
        }

        triangles.insert(triangles.end(), more.begin(), more.end());
    }
};


void d() {
    for (const Point::value_type dx : {0., -2., -1., 1., 1.1}) {
        for (const Point::value_type dy : {0., -2., -1., 1., 1.1}) {
            const Point D{dx, dy};

            const auto A = D + Point{0, 0};
            const auto B = D + Point{M_SQRT2, 0};
            const auto C = D + Point{0, M_SQRT2};

            struct {
                const Point& A;
                const Point& B;
                const Point& C;
            } reorder[] = {
                {A, B, C}, {A, C, B}, {B, C, A}, {B, A, C}, {C, A, B}, {C, B, A},
            };

            //            // should all be the same: circumcenter: {0.707107, 0.707107} circumradius: 1
            //            for (const auto& abc : reorder) {
            //                Triangle T(abc.A, abc.B, abc.C);
            //                std::cout << "Triangle circumcenter: " << (T.circumcenter() - D)
            //                          << " circumradius: " << T.circumradius() << std::endl;
            //            }
        }
    }
}
