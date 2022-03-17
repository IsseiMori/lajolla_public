#include "load_hair.h"
#include <fstream>
#include <iostream>


BezierCurve load_hair(const fs::path &filename, float radius) {
    std::fstream fs(filename.c_str(), std::fstream::in);
    std::string line;
    BezierCurve curve;
    curve.radius = radius;
    while (std::getline(fs, line))
    {
        std::string token;
        std::string delim = " ";
        std::vector<std::string> token_list;
        size_t start = 0;
        size_t end = 0;
        while ((end = line.find(delim, start)) != std::string::npos) {
            token = line.substr(start, end - start);
            // std::cout << token << std::endl;
            token_list.push_back(token);
            start = end + delim.length();
        }

        // read control point and width and ignore everything else
        // control point
        curve.points.push_back(Vector3{std::stof(token_list[10]), std::stof(token_list[11]), std::stof(token_list[12])});
        curve.points.push_back(Vector3{std::stof(token_list[13]), std::stof(token_list[14]), std::stof(token_list[15])});
        curve.points.push_back(Vector3{std::stof(token_list[16]), std::stof(token_list[17]), std::stof(token_list[18])});
        curve.points.push_back(Vector3{std::stof(token_list[19]), std::stof(token_list[20]), std::stof(token_list[21])});
    }
    curve.num = curve.points.size() / 4;
    std::cout << curve.num << "\n";
    return curve;
}

// BezierCurve load_hair(const fs::path &filename, float radius) {
//     std::fstream fs(filename.c_str(), std::fstream::in);
//     BezierCurve curve;
//     curve.radius = radius;
//     // ignore magic number
//     fs.ignore(11);

//     // ignore number
//     fs.ignore(4);
//     float data;
//     int num = 0;
//     Vector3 point{0.0, 0.0, 0.0};
//     while (fs.read((char *)&data, sizeof(float)))
//     {
//         if (data == std::numeric_limits<float>::infinity()) {
//             continue;
//         }
//         point[num] = data;
//         num += 1;
//         if (num == 3) {
//             num = 0;
//             curve.points.push_back(point);
//         }
//     }
//     curve.num = curve.points.size() / 4;
//     std::cout << filename << " loaded: " << curve.num << " curves.\n";
//     return curve;
// }