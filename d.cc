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
#include <cmath>
#include <iostream>


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


struct Triangle {
    Triangle(const Point& A, const Point& B, const Point& C) : A_(A), B_(B), C_(C), circumcenter_(A) {
        auto AB = B - A;
        auto AC = C - A;
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

    const Point& A() const { return A_; }
    const Point& B() const { return B_; }
    const Point& C() const { return C_; }
    const Point& circumcenter() const { return circumcenter_; }
    Point::value_type circumradius() const { return circumradius_; }
    Point::value_type circumradius2() const { return circumradius2_; }

private:
    Point A_;
    Point B_;
    Point C_;
    Point circumcenter_;
    Point::value_type circumradius_;
    Point::value_type circumradius2_;

    friend std::ostream& operator<<(std::ostream& out, const Triangle& T) {
        return out << "{A: " << T.A_ << ", B: " << T.B_ << ", C: " << T.C_ << ", circumcenter: " << T.circumcenter_
                   << ", circumradius: " << T.circumradius_ << '}';
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

            // should all be the same: circumcenter: {0.707107, 0.707107} circumradius: 1
            for (const auto& abc : reorder) {
                Triangle T(abc.A, abc.B, abc.C);
                std::cout << "Triangle circumcenter: " << (T.circumcenter() - D)
                          << " circumradius: " << T.circumradius() << std::endl;
            }
        }
    }
}
