#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define G 39.478// 引力常数（AU^3/(M_sun*yr^2)）
// 适宜轨道距离范围（AU，模拟“宜居带”）
#define MIN_SUIT_DIST 0.8
#define MAX_SUIT_DIST 1.2
// 标准年：行星绕单星公转一周的时间（地球年）
#define STANDARD_YEAR 1.0

// 天体结构体（恒星/行星）
typedef struct {
    double mass;       // 质量（太阳质量M_sun）
    double x, y, z;    // 位置（天文单位AU）
    double vx, vy, vz; // 速度（AU/年）
} CelestialBody;

// 计算两个天体间的引力加速度
void calculateAcceleration(CelestialBody *a, CelestialBody *b, double *ax, double *ay, double *az) {
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    double dz = b->z - a->z;
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    double distCube = dist * dist * dist;
    *ax = G * b->mass * dx / distCube;
    *ay = G * b->mass * dy / distCube;
    *az = G * b->mass * dz / distCube;
}

// 模拟天体运动（欧拉法，单步时间：0.01年）
void simulateMotion(CelestialBody stars[3], CelestialBody *planet, int steps) {
    for (int s = 0; s < steps; s++) {
        double ax = 0, ay = 0, az = 0;
        // 计算三颗恒星对行星的引力加速度
        for (int i = 0; i < 3; i++) {
            double axt, ayt, azt;
            calculateAcceleration(planet, &stars[i], &axt, &ayt, &azt);
            ax += axt;
            ay += ayt;
            az += azt;
        }
        // 更新行星速度（欧拉法）
        double dt = 0.01;
        planet->vx += ax * dt;
        planet->vy += ay * dt;
        planet->vz += az * dt;
        // 更新行星位置
        planet->x += planet->vx * dt;
        planet->y += planet->vy * dt;
        planet->z += planet->vz * dt;
    }
}

// 判断纪元类型：恒纪元/乱纪元
char* judgeEra(CelestialBody stars[3], CelestialBody *planet) {
    // 计算行星与最近恒星的距离
    double minDist = 1e9;
    for (int i = 0; i < 3; i++) {
        double dx = planet->x - stars[i].x;
        double dy = planet->y - stars[i].y;
        double dz = planet->z - stars[i].z;
        double dist = sqrt(dx*dx + dy*dy + dz*dz);
        if (dist < minDist) minDist = dist;
    }
    // 适宜距离内为恒纪元，否则乱纪元
    if (minDist >= MIN_SUIT_DIST && minDist <= MAX_SUIT_DIST) {
        return "恒纪元";
    } else {
        return "乱纪元";
    }
}

// 生成万年历：输出指定年份的纪元类型
void generateCalendar(CelestialBody stars[3], CelestialBody *planet, int years) {
    srand(time(0));
    // 初始化行星轨道（随机位置+速度）
    planet->x = (rand() % 10 - 5) * 0.2;
    planet->y = (rand() % 10 - 5) * 0.2;
    planet->z = (rand() % 10 - 5) * 0.2;
    planet->vx = (rand() % 10 - 5) * 0.1;
    planet->vy = (rand() % 10 - 5) * 0.1;
    planet->vz = (rand() % 10 - 5) * 0.1;

    // 按标准年模拟，输出每年的纪元
    for (int y = 1; y <= years; y++) {
        // 模拟一个标准年的运动（步数=1/0.01=100步）
        simulateMotion(stars, planet, 100);
        char *era = judgeEra(stars, planet);
        printf("第%d年：%s\n", y, era);
    }
}

int main() {
    // 初始化三颗三体恒星（自定义质量和初始位置）
    CelestialBody stars[3] = {
        {1.0, 0, 0, 0, 0, 0, 0},    // 恒星1：质量1M_sun，原点
        {1.0, 2, 0, 0, 0, 0.5, 0},  // 恒星2：质量1M_sun，位置(2,0,0)，速度(0,0.5,0)
        {1.0, -1, 1, 0, 0, -0.3, 0} // 恒星3：质量1M_sun，位置(-1,1,0)，速度(0,-0.3,0)
    };
    // 初始化行星（质量忽略，三体行星）
    CelestialBody planet = {0.001, 1, 0, 0, 0, 2, 0}; // 质量0.001M_sun，初始位置(1,0,0)

    // 生成100年的三体万年历
    printf("三体人万年历（1-100年）：\n");
    generateCalendar(stars, &planet, 100);

    return 0;
}
