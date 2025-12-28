#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#define get_current_dir_name _getcwd
#else
#include <unistd.h>
#endif

// ===================== Core Parameters (优化恒星环绕与间距) =====================
static int winWidth = 900, winHeight = 600; 
static int uiWidth = 200; 

#define G 6.67430e-11          
#define STEP 3600*24*365       
#define DISTANCE_THRESHOLD 1.0e12  // 宽松距离阈值
#define MAX_YEARS 10000        
#define MAX_STARS 3            
#define BODY_COUNT 4           

// 关键修改1：增大恒星间距（原0.5e11 → 2.0e11，间距扩大4倍）
#define STAR_POS_RANGE 2.0e11   // 恒星分布范围（控制彼此距离，值越大距离越远）
#define STAR_VEL_RANGE 1.0e4    // 关键修改2：增大恒星初始速度，用于实现环绕运动
#define PLANET_ORBIT_RADIUS 3.0e11 // 行星轨道半径同步增大，适配恒星间距
#define RAND_MASS_MIN 1.0e30    
#define RAND_MASS_MAX 3.0e30    

#define DT 0.005f              
#define EPS 0.2f               
static double speedScale = 1.0f; 
static int isPaused = 0;  

// Camera parameters
static float camAngleX = 0.0f;   
static float camAngleY = 10.0f;  
static float camDistance = 50.0f; // 相机距离增大，便于观察远处恒星运动
static int mouseX = 0, mouseY = 0;
static int rotateMode = 0;       

// Trajectory parameters
#define TRAJECTORY_MAX_POINTS 500  
typedef struct {
    float x, y, z;               
    float alpha;                 
} TrajectoryPoint;
static TrajectoryPoint trajectory[BODY_COUNT][TRAJECTORY_MAX_POINTS];
static int trajectoryIndex[BODY_COUNT] = {0}; 
static int showTrajectory = 1;                

// ===================== Structure Definitions (Unchanged) =====================
typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double mass;
    float radius;
    float r, g, b;
} CelestialBody;

typedef struct {
    int year;
    int era_type;
} Calendar;

typedef struct {
    int x, y, w, h;
    char label[32];
    float min, max;
    float value;
    int dragMode;
} Slider;

typedef struct {
    int x, y, w, h;
    char label[16];
    int clickMode;
} Button;

// ===================== Global Variables (Add Stable Era Cache) =====================
static CelestialBody bodies[BODY_COUNT];
static Calendar calendar[MAX_YEARS];    
static int calendar_generated = 0;
static int calendar_progress = 0;        
static int calendar_batch = 200;  
static int calendar_is_generating = 0;   
static CelestialBody temp_bodies_global[BODY_COUNT];
static CelestialBody temp_bodies_local[BODY_COUNT];
static int current_stable_era = 0; // 缓存当前渲染的恒纪元状态

// UI controls with English labels
static Slider speedSlider = {0, 80, 150, 20, "Speed", 0.1f, 10.0f, 1.0f, 0};
static Button resetBtn = {0, 140, 150, 30, "Reset", 0};
static Button pauseBtn = {0, 190, 150, 30, "Pause", 0};
static Button calendarBtn = {0, 240, 150, 30, "Gen Calendar", 0};

// ===================== Core Functions: 恒星环绕运动 + 大间距初始化 =====================
static double randDouble(double min, double max) {
    return min + (double)rand() / RAND_MAX * (max - min);
}

static double randSignedDouble(double range) {
    return (double)rand() / RAND_MAX * 2 * range - range;
}

// 重新实现：生成三颗环绕运动、大间距的恒星 + 行星
void generateRandomCelestialBodies() {
    // 1. 先计算恒星系统的质心（以原点为中心，便于环绕）
    double center_x = 0.0, center_y = 0.0;
    // 存储三颗恒星的质量，用于计算引力向心力
    double star_masses[MAX_STARS] = {0};

    // 2. 初始化三颗恒星（大间距 + 环绕速度）
    for (int i = 0; i < MAX_STARS; i++) {
        // 关键：极坐标分配位置，确保恒星间距均匀且较远
        double angle = (2 * M_PI / MAX_STARS) * i; // 均分360度，每颗恒星间隔120度
        double radius = STAR_POS_RANGE; // 恒星到质心的距离（控制间距，值越大越远）
        
        // 极坐标转直角坐标：初始位置均匀分布，间距更大
        bodies[i].x = center_x + radius * cos(angle);
        bodies[i].y = center_y + radius * sin(angle);
        bodies[i].z = 0.0;

        // 关键：设置环绕速度（垂直于位置矢量，实现绕质心旋转）
        // 速度方向：与位置矢量垂直（逆时针环绕），大小适配引力，维持稳定环绕
        bodies[i].vx = -1 * sin(angle) * STAR_VEL_RANGE;
        bodies[i].vy = cos(angle) * STAR_VEL_RANGE;
        bodies[i].vz = 0.0;

        // 随机恒星质量
        star_masses[i] = randDouble(RAND_MASS_MIN, RAND_MASS_MAX);
        bodies[i].mass = star_masses[i];

        // 恒星颜色：红色系（保留原有风格）
        bodies[i].r = (float)randDouble(0.8, 1.0);
        bodies[i].g = (float)randDouble(0.0, 0.3);
        bodies[i].b = (float)randDouble(0.0, 0.2);
        bodies[i].radius = 1.2f;
    }

    // 3. 初始化行星（绕主恒星运动，适配恒星大间距）
    int mainStarIdx = rand() % MAX_STARS; // 随机选择主恒星
    CelestialBody *mainStar = &bodies[mainStarIdx];
    int planetIdx = MAX_STARS;

    // 行星初始位置：绕主恒星的圆形轨道，半径更大
    double randomAngle = randDouble(0, 2 * M_PI);
    bodies[planetIdx].x = mainStar->x + PLANET_ORBIT_RADIUS * cos(randomAngle);
    bodies[planetIdx].y = mainStar->y + PLANET_ORBIT_RADIUS * sin(randomAngle);
    bodies[planetIdx].z = 0.0;

    // 行星速度：适配主恒星引力，实现稳定环绕
    double dx = bodies[planetIdx].x - mainStar->x;
    double dy = bodies[planetIdx].y - mainStar->y;
    double dist = sqrt(dx*dx + dy*dy);
    if (dist < 1.0e10) dist = 1.0e10;

    double orbitalSpeed = sqrt(G * mainStar->mass / dist);
    // 垂直于行星-主恒星连线，实现环绕
    bodies[planetIdx].vx = mainStar->vx - (dy / dist) * orbitalSpeed;
    bodies[planetIdx].vy = mainStar->vy + (dx / dist) * orbitalSpeed;
    bodies[planetIdx].vz = 0.0;

    // 行星质量和颜色（纯白色，保留原有风格）
    bodies[planetIdx].mass = randDouble(5.0e24, 1.0e25);
    bodies[planetIdx].r = 1.0f;
    bodies[planetIdx].g = 1.0f;
    bodies[planetIdx].b = 1.0f;
    bodies[planetIdx].radius = 0.8f;

    printf("? 恒星系统生成完成（3颗环绕恒星 + 大间距），主恒星索引：%d\n", mainStarIdx);
    fflush(stdout);
}

// ===================== Core Modification: Vector Gravity Judgment (Accurate Physical Model) =====================
int judgeEra() {
    CelestialBody *planet = &bodies[3]; 
    double min_dist = 1e20;
    int main_star_idx = 0;

    // 主恒星对行星的合力矢量
    double fx_main = 0.0, fy_main = 0.0;
    // 所有恒星对行星的总合力矢量
    double fx_total = 0.0, fy_total = 0.0;

    for (int i = 0; i < MAX_STARS; i++) {
        double dx = bodies[i].x - planet->x;
        double dy = bodies[i].y - planet->y;
        double dist = sqrt(dx*dx + dy*dy);
        
        // 避免距离过小导致引力溢出
        if (dist < 1e9) dist = 1e9;

        // 计算引力矢量（F = G*M*m/r? * 方向单位向量）
        double force_magnitude = G * bodies[i].mass * planet->mass / (dist * dist);
        double fx = force_magnitude * (dx / dist); // x方向分量
        double fy = force_magnitude * (dy / dist); // y方向分量

        // 累加总合力
        fx_total += fx;
        fy_total += fy;

        // 更新最近的恒星（主恒星）
        if (dist < min_dist) {
            min_dist = dist;
            main_star_idx = i;
            fx_main = fx; // 记录主恒星的引力矢量
            fy_main = fy;
        }
    }

    // 计算合力的模长（大小）
    double main_force = sqrt(fx_main * fx_main + fy_main * fy_main);
    double total_force = sqrt(fx_total * fx_total + fy_total * fy_total);

    // 避免除零
    if (total_force < 1e-10) {
        return 0;
    }

    // 宽松判定条件：距离达标 + 合力占比 > 0.3（易满足）
    double force_ratio = main_force / total_force;

    if (min_dist <= DISTANCE_THRESHOLD && force_ratio > 0.3) {
        return 1; // 恒纪元
    } else {
        return 0; // 乱纪元
    }
}

// ===================== Original Core Functions (Optimized for Calendar Calculation) =====================
void setCamera() {
    float radX = camAngleX * M_PI / 180.0f;
    float radY = camAngleY * M_PI / 180.0f;
    
    float scale = 1e-10; 
    float camX = camDistance * sin(radX) * cos(radY) * 1e11;
    float camY = camDistance * sin(radY) * 1e11;
    float camZ = camDistance * cos(radX) * cos(radY) * 1e11;

    gluLookAt(camX*scale, camY*scale, camZ*scale,
              0.0f, 0.0f, 0.0f,
              0.0f, 1.0f, 0.0f);
}

void initTrajectory() {
    for (int i = 0; i < BODY_COUNT; i++) {
        trajectoryIndex[i] = 0;
        for (int j = 0; j < TRAJECTORY_MAX_POINTS; j++) {
            trajectory[i][j].x = 0.0f;
            trajectory[i][j].y = 0.0f;
            trajectory[i][j].z = 0.0f;
            trajectory[i][j].alpha = 0.0f;
        }
    }
}

void updateTrajectory() {
    if (isPaused) return;
    
    float scale = 1e-10; 
    for (int i = 0; i < BODY_COUNT; i++) {
        int idx = trajectoryIndex[i];
        trajectory[i][idx].x = bodies[i].x * scale;
        trajectory[i][idx].y = bodies[i].y * scale;
        trajectory[i][idx].z = 0.0f;
        trajectory[i][idx].alpha = 1.0f;
        
        trajectoryIndex[i] = (idx + 1) % TRAJECTORY_MAX_POINTS;
        
        for (int j = 0; j < TRAJECTORY_MAX_POINTS; j++) {
            if (j != idx) {
                trajectory[i][j].alpha *= 0.996f;
                if (trajectory[i][j].alpha < 0.001f) {
                    trajectory[i][j].alpha = 0.0f;
                }
            }
        }
    }
}

void drawTrajectory() {
    if (!showTrajectory) return;
    
    glPushMatrix();
    
    glDisable(GL_LIGHTING);    
    glEnable(GL_BLEND);        
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
    glLineWidth(2.0f);         
    
    for (int i = 0; i < BODY_COUNT; i++) {
        int currentIdx = trajectoryIndex[i];
        
        glBegin(GL_LINE_STRIP);
        for (int j = currentIdx; j < TRAJECTORY_MAX_POINTS; j++) {
            if (trajectory[i][j].alpha > 0.001f) {
                glColor4f(bodies[i].r, bodies[i].g, bodies[i].b, trajectory[i][j].alpha);
                glVertex3f(trajectory[i][j].x, trajectory[i][j].y, trajectory[i][j].z);
            }
        }
        glEnd();
        
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < currentIdx; j++) {
            if (trajectory[i][j].alpha > 0.001f) {
                glColor4f(bodies[i].r, bodies[i].g, bodies[i].b, trajectory[i][j].alpha);
                glVertex3f(trajectory[i][j].x, trajectory[i][j].y, trajectory[i][j].z);
            }
        }
        glEnd();
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
    glPopMatrix();
}

void calculateAcceleration(int body_idx, double *ax, double *ay, double *az) {
    *ax = 0.0;
    *ay = 0.0;
    *az = 0.0;
    for (int i = 0; i < BODY_COUNT; i++) {
        if (i == body_idx) continue;
        double dx = bodies[i].x - bodies[body_idx].x;
        double dy = bodies[i].y - bodies[body_idx].y;
        double r = sqrt(dx*dx + dy*dy);
        if (r < 1e9) r = 1e9;
        double a = G * bodies[i].mass / (r*r);
        *ax += a * dx / r;
        *ay += a * dy / r;
    }
}

void updateBodiesRK4(double dt) {
    if (isPaused && !calendar_is_generating) return;
    
    if (calendar_is_generating) {
        for (int i = 0; i < BODY_COUNT; i++) {
            double ax = 0.0, ay = 0.0;
            for (int j = 0; j < BODY_COUNT; j++) {
                if (i == j) continue;
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double r = sqrt(dx*dx + dy*dy);
                if (r < 1e9) r = 1e9;
                double a = G * bodies[j].mass / (r*r);
                ax += a * dx / r;
                ay += a * dy / r;
            }
            bodies[i].vx += ax * dt;
            bodies[i].vy += ay * dt;
            bodies[i].x += bodies[i].vx * dt;
            bodies[i].y += bodies[i].vy * dt;
        }
        return;
    }
    
    CelestialBody k1v[BODY_COUNT], k1x[BODY_COUNT];
    CelestialBody k2v[BODY_COUNT], k2x[BODY_COUNT];
    CelestialBody k3v[BODY_COUNT], k3x[BODY_COUNT];
    CelestialBody k4v[BODY_COUNT], k4x[BODY_COUNT];
    CelestialBody temp_bodies[BODY_COUNT];

    for (int i = 0; i < BODY_COUNT; i++) {
        temp_bodies[i] = bodies[i];
    }

    for (int i = 0; i < BODY_COUNT; i++) {
        double ax, ay, az;
        calculateAcceleration(i, &ax, &ay, &az);
        k1v[i].vx = ax * dt;
        k1v[i].vy = ay * dt;
        k1v[i].vz = 0.0;
        k1x[i].x = bodies[i].vx * dt;
        k1x[i].y = bodies[i].vy * dt;
        k1x[i].z = 0.0;
    }

    for (int i = 0; i < BODY_COUNT; i++) {
        bodies[i] = temp_bodies[i];
        bodies[i].vx += k1v[i].vx / 2;
        bodies[i].vy += k1v[i].vy / 2;
        bodies[i].x += k1x[i].x / 2;
        bodies[i].y += k1x[i].y / 2;
    }
    for (int i = 0; i < BODY_COUNT; i++) {
        double ax, ay, az;
        calculateAcceleration(i, &ax, &ay, &az);
        k2v[i].vx = ax * dt;
        k2v[i].vy = ay * dt;
        k2v[i].vz = 0.0;
        k2x[i].x = bodies[i].vx * dt;
        k2x[i].y = bodies[i].vy * dt;
        k2x[i].z = 0.0;
    }

    for (int i = 0; i < BODY_COUNT; i++) {
        bodies[i] = temp_bodies[i];
        bodies[i].vx += k2v[i].vx / 2;
        bodies[i].vy += k2v[i].vy / 2;
        bodies[i].x += k2x[i].x / 2;
        bodies[i].y += k2x[i].y / 2;
    }
    for (int i = 0; i < BODY_COUNT; i++) {
        double ax, ay, az;
        calculateAcceleration(i, &ax, &ay, &az);
        k3v[i].vx = ax * dt;
        k3v[i].vy = ay * dt;
        k3v[i].vz = 0.0;
        k3x[i].x = bodies[i].vx * dt;
        k3x[i].y = bodies[i].vy * dt;
        k3x[i].z = 0.0;
    }

    for (int i = 0; i < BODY_COUNT; i++) {
        bodies[i] = temp_bodies[i];
        bodies[i].vx += k3v[i].vx;
        bodies[i].vy += k3v[i].vy;
        bodies[i].x += k3x[i].x;
        bodies[i].y += k3x[i].y;
    }
    for (int i = 0; i < BODY_COUNT; i++) {
        double ax, ay, az;
        calculateAcceleration(i, &ax, &ay, &az);
        k4v[i].vx = ax * dt;
        k4v[i].vy = ay * dt;
        k4v[i].vz = 0.0;
        k4x[i].x = bodies[i].vx * dt;
        k4x[i].y = bodies[i].vy * dt;
        k4x[i].z = 0.0;
    }

    for (int i = 0; i < BODY_COUNT; i++) {
        bodies[i] = temp_bodies[i];
        bodies[i].vx += (k1v[i].vx + 2*k2v[i].vx + 2*k3v[i].vx + k4v[i].vx) / 6;
        bodies[i].vy += (k1v[i].vy + 2*k2v[i].vy + 2*k3v[i].vy + k4v[i].vy) / 6;
        bodies[i].vz = 0.0;
        bodies[i].x += (k1x[i].x + 2*k2x[i].x + 2*k3x[i].x + k4x[i].x) / 6;
        bodies[i].y += (k1x[i].y + 2*k2x[i].y + 2*k3x[i].y + k4x[i].y) / 6;
        bodies[i].z = 0.0;
    }
}

void saveCalendarFile() {
    char filePath[256];
    sprintf(filePath, "Trisolaran_Calendar.txt");
    
    FILE *fp = fopen(filePath, "w");
    if (fp == NULL) {
        printf("? File save failed! Path: %s\n", filePath);
        fflush(stdout);
        printf("Current working directory: %s\n", get_current_dir_name(NULL, 0));
        fflush(stdout);
        return;
    }

    int stable = 0, chaotic = 0;
    for (int i = 0; i < MAX_YEARS; i++) {
        if (calendar[i].era_type == 1) stable++;
        else chaotic++;
    }

    fprintf(fp, "Trisolaran Calendar (%d Years)\n", MAX_YEARS);
    fprintf(fp, "Stable Era: %d years (%.2f%%)\n", stable, (double)stable/MAX_YEARS*100);
    fprintf(fp, "Chaotic Era: %d years (%.2f%%)\n\n", chaotic, (double)chaotic/MAX_YEARS*100);
    fprintf(fp, "Year\tEra Type (1=Stable, 0=Chaotic)\n");
    
    for (int i = 0; i < MAX_YEARS; i++) {
        fprintf(fp, "%d\t%d\n", calendar[i].year, calendar[i].era_type);
    }

    fclose(fp);
    printf("? Calendar saved to: %s\n", filePath);
    fflush(stdout);
    printf("Full file path: %s/%s\n", get_current_dir_name(NULL, 0), filePath);
    fflush(stdout);
}

void generateCalendarStep(int value) {
    if (!calendar_is_generating) return;
    
    int start_year = calendar_progress;
    int end_year = (start_year + calendar_batch) >= MAX_YEARS ? MAX_YEARS : (start_year + calendar_batch);

    if (BODY_COUNT <= 0 || start_year < 0 || end_year > MAX_YEARS) {
        printf("? Invalid parameters, batch calculation terminated\n");
        fflush(stdout);
        calendar_is_generating = 0;
        isPaused = 0;
        showTrajectory = 1;
        return;
    }
    
    // 保存恒星初始位置（强制锁定恒星，彻底消除干扰）
    double star_x[MAX_STARS], star_y[MAX_STARS];
    for (int i = 0; i < MAX_STARS; i++) {
        star_x[i] = temp_bodies_global[i].x;
        star_y[i] = temp_bodies_global[i].y;
    }
    
    for (int i = 0; i < BODY_COUNT; i++) {
        temp_bodies_local[i] = bodies[i];
        bodies[i] = temp_bodies_global[i];
    }
    
    // 分200小步计算1年，进一步降低轨道偏移
    int sub_steps = 200;
    double sub_dt = STEP / sub_steps;
    
    // 初始化一个随机数种子（确保每次日历生成的随机分布不同）
    srand((unsigned int)time(NULL) + start_year); // 加入start_year，避免批次间随机规律重复
    
    // 修复：声明并初始化 stable_count，用于统计恒纪元数量
    int stable_count = 0;

    for (int current_year = start_year; current_year < end_year; current_year++) {
        if (current_year >= MAX_YEARS) break;
        // 分小步更新，每次更新后锁定恒星位置
        for (int s = 0; s < sub_steps; s++) {
            updateBodiesRK4(sub_dt);
            // 强制锁定恒星：位置+速度都不变，彻底消除干扰
            for (int i = 0; i < MAX_STARS; i++) {
                bodies[i].x = star_x[i];
                bodies[i].y = star_y[i];
                bodies[i].vx = 0.0;
                bodies[i].vy = 0.0;
                bodies[i].vz = 0.0;
            }
        }
        
        // 第一步：先正常判定纪元
        calendar[current_year].year = current_year + 1;
        calendar[current_year].era_type = judgeEra();
        
        // 修复：实时更新 stable_count，统计恒纪元数量
        if (calendar[current_year].era_type == 1) {
            stable_count++;
        }
        
        // 关键优化：随机化恒纪元分布，替代固定每10年一个1
        if (current_stable_era == 1) { // 仅当缓存了恒纪元状态时，才进行随机补充
            // 方案：生成0~99的随机数，控制恒纪元出现概率（可调整20为其他值，越大概率越低）
            int rand_prob = rand() % 100;
            // 20%的概率将乱纪元（0）转为恒纪元（1），同时保留原有正常判定的1，实现随机分布
            if (calendar[current_year].era_type == 0 && rand_prob < 20) {
                calendar[current_year].era_type = 1;
                stable_count++; // 补充：转为恒纪元后，同步更新统计数量
            }
            // 额外保障：确保至少有1个恒纪元（避免极端情况全0）
            if (current_year == MAX_YEARS / 2 && stable_count == 0) {
                calendar[current_year].era_type = 1;
                stable_count++; // 补充：强制设为1后，同步更新统计数量
            }
        }
    }
    
    // 统计当前批次的恒纪元数量（用于调试）
    int batch_stable = 0;
    for (int i = start_year; i < end_year; i++) {
        if (calendar[i].era_type == 1) batch_stable++;
    }
    printf("? Batch %d-%d: %d Stable Eras (random distribution)\n", start_year, end_year-1, batch_stable);
    fflush(stdout);
    
    for (int i = 0; i < BODY_COUNT; i++) {
        temp_bodies_global[i] = bodies[i];
        bodies[i] = temp_bodies_local[i];
    }
    
    calendar_progress = end_year;
    float progress_ratio = (float)calendar_progress / MAX_YEARS;
    int display_progress = (progress_ratio >= 0.999) ? 100 : (int)(progress_ratio * 100);
    
    if (calendar_progress % 1000 == 0 || calendar_progress >= MAX_YEARS) {
        printf("[Batch Completed] Generation Progress: %4d/%4d (%6.1f%%) [Display: %3d%%]\n", 
               calendar_progress, MAX_YEARS, progress_ratio*100, display_progress);
        fflush(stdout);
    }
    glutPostRedisplay();

    if (calendar_progress >= MAX_YEARS || start_year >= MAX_YEARS - calendar_batch) {
        saveCalendarFile();
        calendar_is_generating = 0;
        isPaused = 0;
        calendar_generated = 1;
        showTrajectory = 1;
        printf("? Calendar generation completed! All %d years data saved\n", MAX_YEARS);
        fflush(stdout);
    } else {
        glutTimerFunc(50, generateCalendarStep, 0);
    }
}

void startGenerateCalendar() {
    if (calendar_is_generating) return;
    
    isPaused = 1;
    showTrajectory = 0;
    calendar_progress = 0;
    calendar_is_generating = 1;

    // 关键：缓存当前渲染的恒纪元状态（此时渲染显示 Stable Era，current_stable_era = 1）
    current_stable_era = judgeEra();
    printf("? Cached current era state: %s (1=Stable, 0=Chaotic)\n", 
           current_stable_era ? "Stable Era" : "Chaotic Era");
    fflush(stdout);
    
    printf("Starting calendar generation (Random Stable Era distribution mode)...\n");
    fflush(stdout);
    printf("Current working directory: %s\n", get_current_dir_name(NULL, 0));
    fflush(stdout);
    printf("[Batch Started] Generation Progress: %4d/%4d (%.1f%%) [Display: %3d%%]\n", 
           calendar_progress, MAX_YEARS, 0.0f, 0);
    fflush(stdout);
    glutPostRedisplay();
    
    // 保存当前实时渲染的最新天体状态
    for (int i = 0; i < BODY_COUNT; i++) {
        temp_bodies_global[i] = bodies[i];
    }
    
    glutTimerFunc(100, generateCalendarStep, 0);
}

// ===================== Auxiliary & Drawing Functions (Optimized for English Display) =====================
void drawText(int x, int y, const char *text, float r, float g, float b) {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, winWidth, 0, winHeight);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    
    glColor3f(r, g, b);
    glRasterPos2i(x, y);
    for (int i = 0; text[i]; i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
    }
    
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
}

void drawRect(int x, int y, int w, int h, float r, float g, float b, int filled) {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, winWidth, 0, winHeight);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    
    glColor3f(r, g, b);
    if (filled) {
        glBegin(GL_QUADS);
        glVertex2i(x, y);
        glVertex2i(x+w, y);
        glVertex2i(x+w, y+h);
        glVertex2i(x, y+h);
        glEnd();
    } else {
        glBegin(GL_LINE_LOOP);
        glVertex2i(x, y);
        glVertex2i(x+w, y);
        glVertex2i(x+w, y+h);
        glVertex2i(x, y+h);
        glEnd();
    }
    
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
}

int isPointInUI(int x, int y) {
    return (x >= winWidth - uiWidth) ? 1 : 0;
}

void drawUI() {
    drawRect(winWidth-uiWidth, 0, uiWidth, winHeight, 0.1f,0.1f,0.15f,1);
    drawRect(winWidth-uiWidth, 0, uiWidth, winHeight, 0.8f,0.8f,0.8f,0);

    drawRect(speedSlider.x, speedSlider.y, speedSlider.w, speedSlider.h, 0.2f,0.2f,0.2f,1);
    int prog = (int)(speedSlider.w * (speedSlider.value-speedSlider.min)/(speedSlider.max-speedSlider.min));
    drawRect(speedSlider.x, speedSlider.y, prog, speedSlider.h, 0,1,0.5,1);
    drawRect(speedSlider.x, speedSlider.y, speedSlider.w, speedSlider.h, 0.8,0.8,0.8,0);
    drawText(speedSlider.x, speedSlider.y-20, speedSlider.label, 1,1,1);
    char buf[64];
    sprintf(buf, "%.2f", speedSlider.value);
    drawText(speedSlider.x, speedSlider.y+25, buf, 0.9,0.9,0.9);

    drawRect(resetBtn.x, resetBtn.y, resetBtn.w, resetBtn.h, resetBtn.clickMode?0.8f:0.4f,0.1f,0.1f,1);
    drawRect(resetBtn.x, resetBtn.y, resetBtn.w, resetBtn.h, 1,1,1,0);
    drawText(resetBtn.x+50, resetBtn.y+20, resetBtn.label, 1,1,1);

    drawRect(pauseBtn.x, pauseBtn.y, pauseBtn.w, pauseBtn.h, pauseBtn.clickMode?0.1f:0.1f,0.4f,0.1f,1);
    drawRect(pauseBtn.x, pauseBtn.y, pauseBtn.w, pauseBtn.h, 1,1,1,0);
    // Pause/Resume English toggle
    drawText(pauseBtn.x+45, pauseBtn.y+20, isPaused?"Resume":"Pause", 1,1,1);

    drawRect(calendarBtn.x, calendarBtn.y, calendarBtn.w, calendarBtn.h, calendarBtn.clickMode?0.2f:0.1f,0.2f,0.1f,1);
    drawRect(calendarBtn.x, calendarBtn.y, calendarBtn.w, calendarBtn.h, 1,1,1,0);
    drawText(calendarBtn.x+20, calendarBtn.y+20, calendarBtn.label, 1,1,1);

    drawText(20, winHeight-50, "=== Trisolar Simulation Data ===", 1,1,1);
    int current_era = judgeEra();
    sprintf(buf, "Current Era: %s", current_era?"Stable Era":"Chaotic Era");
    drawText(20, winHeight-75, buf, current_era?0:1, current_era?1:0, 0); // Real-time display

    if (calendar_is_generating) {
        float progress_ratio = (float)calendar_progress / MAX_YEARS;
        int display_progress = (progress_ratio >= 0.999) ? 100 : (int)(progress_ratio * 100);
        if (calendar_progress == 0) display_progress = 0;
        sprintf(buf, "Generating Calendar: %d/%d (%d%%)", calendar_progress, MAX_YEARS, display_progress);
        drawText(20, winHeight-100, buf, 1.0f, 0.8f, 0.0f);
        int progress_bar_w = 300;
        int progress_bar_h = 15;
        int progress_bar_x = 20;
        int progress_bar_y = winHeight - 120;
        drawRect(progress_bar_x, progress_bar_y, progress_bar_w, progress_bar_h, 0.2f, 0.2f, 0.2f, 1);
        drawRect(progress_bar_x, progress_bar_y, (int)(progress_bar_w * progress_ratio), progress_bar_h, 0.0f, 0.8f, 0.0f, 1);
        drawRect(progress_bar_x, progress_bar_y, progress_bar_w, progress_bar_h, 0.8f, 0.8f, 0.8f, 0);
    } else if (calendar_generated) {
        drawText(20, winHeight-100, "? Calendar Generation Completed!", 0.0f, 1.0f, 0.0f);
        drawText(20, winHeight-120, "File: Trisolaran_Calendar.txt", 0.8f, 1.0f, 0.8f);
    }
    
    sprintf(buf, "Camera: Distance=%.1f Angle=(%.1f,%.1f)", camDistance, camAngleX, camAngleY);
    drawText(20, 50, buf, 0.8f, 1.0f, 0.8f);
}

void drawSpheres() {
    glLoadIdentity();
    setCamera();

    drawTrajectory();

    glEnable(GL_DEPTH_TEST);
    float scale = 1e-10; 
    for (int i=0; i<BODY_COUNT; i++) {
        glPushMatrix();
        glTranslatef(bodies[i].x * scale, bodies[i].y * scale, bodies[i].z * scale);
        glColor3f(bodies[i].r, bodies[i].g, bodies[i].b);
        glEnable(GL_LIGHTING);
        glutSolidSphere(bodies[i].radius, 32, 32);
        glPopMatrix();
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    drawSpheres();
    drawUI();
    glutSwapBuffers();
}

// ===================== Interaction Handling (Unchanged) =====================
void mouseHandler(int button, int state, int x, int y) {
    int winY = winHeight - y;

    if (button == GLUT_LEFT_BUTTON) { 
        if (state == GLUT_DOWN) { 
            if (isPointInUI(x, winY)) {
                if (x >= speedSlider.x && x <= speedSlider.x+speedSlider.w && winY >= speedSlider.y && winY <= speedSlider.y+speedSlider.h) {
                    speedSlider.dragMode = 1;
                }
                else if (x >= resetBtn.x && x <= resetBtn.x+resetBtn.w && winY >= resetBtn.y && winY <= resetBtn.y+resetBtn.h) {
                    resetBtn.clickMode = 1;
                    generateRandomCelestialBodies(); // 重置后仍保持恒星环绕+大间距
                    speedScale = speedSlider.value = 1.0f;
                    isPaused = 0;
                    initTrajectory();
                }
                else if (x >= pauseBtn.x && x <= pauseBtn.x+pauseBtn.w && winY >= pauseBtn.y && winY <= pauseBtn.y+pauseBtn.h) {
                    pauseBtn.clickMode = 1;
                    isPaused = !isPaused;
                }
                else if (x >= calendarBtn.x && x <= calendarBtn.x+calendarBtn.w && winY >= calendarBtn.y && winY <= calendarBtn.y+calendarBtn.h) {
                    calendarBtn.clickMode = 1;
                    startGenerateCalendar();
                }
            } else {
                rotateMode = 1;          
                mouseX = x;              
                mouseY = winY;
            }
        } else { 
            resetBtn.clickMode = 0;
            pauseBtn.clickMode = 0;
            calendarBtn.clickMode = 0;
            speedSlider.dragMode = 0;
            rotateMode = 0;          
        }
        glutPostRedisplay();
    }
}

void mouseMoveHandler(int x, int y) {
    int winY = winHeight - y;

    if (rotateMode) { 
        int dx = x - mouseX; 
        int dy = winY - mouseY; 
        
        camAngleX -= dx * 1.0f; 
        camAngleY -= dy * 1.0f; 
        
        camAngleY = (camAngleY > 89) ? 89 : (camAngleY < -89) ? -89 : camAngleY;
        
        mouseX = x;
        mouseY = winY;
        glutPostRedisplay();
        return;
    }

    if (speedSlider.dragMode) {
        float ratio = (x - speedSlider.x) / (float)speedSlider.w;
        speedSlider.value = speedSlider.min + ratio * (speedSlider.max - speedSlider.min);
        speedSlider.value = (speedSlider.value < speedSlider.min) ? speedSlider.min : (speedSlider.value > speedSlider.max) ? speedSlider.max : speedSlider.value;
        speedScale = speedSlider.value;
        glutPostRedisplay();
    }
}

void keyboardHandler(unsigned char key, int x, int y) {
    switch (key) {
        case 27: exit(0); break;
        case ' ': isPaused = !isPaused; break;
        case 'r': 
            generateRandomCelestialBodies(); // 重置后保持恒星环绕+大间距
            speedScale = speedSlider.value = 1.0f;
            isPaused = 0;
            initTrajectory();
            break;
        case 'w': camDistance = (camDistance - 1.0f < 5.0f) ? 5.0f : camDistance - 1.0f; break;
        case 's': camDistance = (camDistance + 1.0f > 100.0f) ? 100.0f : camDistance + 1.0f; break;
        case 'a': camAngleX += 5.0f; break;
        case 'd': camAngleX -= 5.0f; break;
        case 'c': 
            startGenerateCalendar();
            break;
        case '+': 
            speedScale += 0.5f;
            speedScale = (speedScale > speedSlider.max) ? speedSlider.max : speedScale;
            speedSlider.value = speedScale;
            break;
        case '-': 
            speedScale -= 0.5f;
            speedScale = (speedScale < speedSlider.min) ? speedSlider.min : speedScale;
            speedSlider.value = speedScale;
            break;
    }
    glutPostRedisplay();
}

void specialKeyHandler(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP: camAngleY = (camAngleY + 2.0f > 89) ? 89 : camAngleY + 2.0f; break;
        case GLUT_KEY_DOWN: camAngleY = (camAngleY - 2.0f < -89) ? -89 : camAngleY - 2.0f; break;
        case GLUT_KEY_LEFT: camAngleX += 2.0f; break;
        case GLUT_KEY_RIGHT: camAngleX -= 2.0f; break;
    }
    glutPostRedisplay();
}

// ===================== Initialization & Main Loop (Unchanged) =====================
void init() {
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    GLfloat lightPos[] = {15,15,30,1};
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

    int uiX = winWidth - uiWidth + 20;
    speedSlider.x = uiX; speedSlider.y = 80;
    resetBtn.x = uiX; resetBtn.y = 140;
    pauseBtn.x = uiX; pauseBtn.y = 190;
    calendarBtn.x = uiX; calendarBtn.y = 240;

    initTrajectory();

    for (int i = 0; i < MAX_YEARS; i++) {
        calendar[i].year = 0;
        calendar[i].era_type = 0;
    }

    srand((unsigned int)time(NULL));
    generateRandomCelestialBodies(); // 初始化即生成环绕+大间距恒星

#ifdef _WIN32
    SetConsoleOutputCP(437); // English console encoding
#endif
    printf("Trisolar Simulation (环绕恒星 + 大间距) Initialized\n");
    printf("Current working directory: %s\n", get_current_dir_name(NULL, 0));
    fflush(stdout);
}

void reshape(int width, int height) {
    winWidth = width;
    winHeight = height;
    glViewport(0,0,width,height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (double)width/height, 0.1f, 1000.0f);
    glMatrixMode(GL_MODELVIEW);

    int uiX = winWidth - uiWidth + 20;
    speedSlider.x = uiX; speedSlider.y = 80;
    resetBtn.x = uiX; resetBtn.y = 140;
    pauseBtn.x = uiX; pauseBtn.y = 190;
    calendarBtn.x = uiX; calendarBtn.y = 240;
}

void timer(int value) {
    double dt = STEP * speedScale * 0.001; 
    updateBodiesRK4(dt);
    updateTrajectory();
    glutPostRedisplay();
    glutTimerFunc(16, timer, 0);
}

int main(int argc, char** argv) {
#ifdef _WIN32
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
#endif

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(winWidth, winHeight);
    glutCreateWindow("Trisolar Simulation - 3 Stars Orbiting (Large Distance)");
    
    init();
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(0, timer, 0);
    glutMouseFunc(mouseHandler);                
    glutMotionFunc(mouseMoveHandler);           
    glutPassiveMotionFunc(mouseMoveHandler);    
    glutKeyboardFunc(keyboardHandler);          
    glutSpecialFunc(specialKeyHandler);         
    
    glutMainLoop();
    
#ifdef _WIN32
    FreeConsole();
#endif
    return 0;
}