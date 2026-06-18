#include <boost/multi/array.hpp>

#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>

#include <algorithm>
#include <limits>
#include <cmath>

namespace multi = boost::multi;
namespace gil   = boost::gil;

int main() {

    constexpr int width  = 120;
    constexpr int height = 100;

    multi::array<std::complex<double>, 2> a = [](auto x, auto y) {
        auto z = std::complex<double>(x/120.0, y/100.0);
        return z*z;
    } ^ multi::extents_t<2>({height, width});

    // auto [min_it, max_it] = std::minmax_element(a.elements().begin(), a.elements().end());

    gil::rgb8_image_t img(width, height);
    auto dst = view(img);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            constexpr auto thresh = 0.1;
            constexpr auto f = [](auto z) { return z*z; };

            // atan2 returns a range from -pi to pi, so we need to add pi, but this offsets
            // the value by 180 degrees, so we also imput (-y, -x) for another 180 degrees
            // to invert rotation
            constexpr auto angle = [](auto x, auto y) { return (M_PI + atan2(-y,-x)) / (2*M_PI); };

            constexpr auto r = [](auto x, auto y){ return sqrt(x*x + y*y);};

            // complex phase and magnitude
            constexpr auto theta = [](auto x, auto y) {return atan2(y,x);};
            constexpr auto z = [r, theta](auto x, auto y) {return r(x,y)*exp(theta(x,y)*std::complex<double>(0, 1)); };

            // imaginary and real output functions
            constexpr auto real_f = [f](auto z) { return real(f(z)); };
            constexpr auto imaginary_f = [f](auto z) {return imag(f(z));};

            // magnitude contours
            auto magnitude_shading = [=](auto x, auto y) { return 0.5 + 0.5*(abs(f(z(x,y)))-floor(abs(f(z(x,y))))); };

            // gridlines
            auto gridlines = [=](auto x, auto y) { return
                pow(abs(sin(real_f     (z(x,y))*M_PI)), thresh)* 
                pow(abs(sin(imaginary_f(z(x,y))*M_PI)), thresh);
            };
            // // Normalize to [0,1]
            // float v = (a[y][x] - minv) * scale;
            // v = std::clamp(v, 0.0f, 1.0f);

            float hue = angle(real_f(z(x,y)), imaginary_f(z(x,y)));
            float saturation = magnitude_shading(x,y);
            float value = gridlines(x,y);
            // HSV encoding
            // float hue        = 360.0f * std::atan2(real(a[y][x]), imag(a[y][x])) / M_PI*4;
            // float saturation = 0.5 + 0.5*std::log(std::abs(a[y][x]));  // 1.0f;
            // float value      = 1.0f;

            // Convert HSV → RGB (Boost.GIL)
            gil::rgb32f_pixel_t rgb_f;
            gil::hsv32f_pixel_t hsv_f{hue, saturation, value};

            gil::color_convert(hsv_f, rgb_f);

            // Convert to 8-bit RGB
            dst(x, y) = gil::rgb8_pixel_t{
                static_cast<unsigned char>(std::lround(rgb_f[0] * 255.0f)),
                static_cast<unsigned char>(std::lround(rgb_f[1] * 255.0f)),
                static_cast<unsigned char>(std::lround(rgb_f[2] * 255.0f))
            };
        }
    }

    // 5) Write PNG
    gil::write_view("output.png", dst, gil::png_tag{});
}
