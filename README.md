#include <math.h>
#include <stdio.h>

/*
 * 作业说明：
 * 已知大地坐标（经度、纬度、高程），完成：
 * 1) 大地坐标 -> 空间直角坐标(ECEF)
 * 2) 空间直角坐标 -> 大地坐标（迭代）
 * 3) 在给定纬度处，1" 对应的距离（纬度方向与经度方向）
 * 4) 若直角坐标精度为 1 mm，则角秒应保留几位小数
 */

/* 角度和弧度转换 */
#define DEG2RAD(x) ((x) * M_PI / 180.0)
#define RAD2DEG(x) ((x) * 180.0 / M_PI)

/* WGS-84/CGCS2000 常用椭球参数（两者 a,f 一致） */
static const double A = 6378137.0;                  /* 长半轴 (m) */
static const double F = 1.0 / 298.257223563;        /* 扁率 */
static const double E2 = F * (2.0 - F);             /* 第一偏心率平方 */

/* DMS(度分秒) -> 十进制度 */
double dms_to_deg(double d, double m, double s) {
    return d + m / 60.0 + s / 3600.0;
}

/* 十进制度 -> DMS(度分秒) */
void deg_to_dms(double deg, int *d, int *m, double *s) {
    double abs_deg = fabs(deg);
    *d = (int)abs_deg;
    double rem_min = (abs_deg - *d) * 60.0;
    *m = (int)rem_min;
    *s = (rem_min - *m) * 60.0;

    /* 保留符号在度上 */
    if (deg < 0) {
        *d = -*d;
    }
}

/* 大地坐标(经纬高) -> ECEF(XYZ) */
void geodetic_to_ecef(double lon_deg, double lat_deg, double h,
                      double *X, double *Y, double *Z) {
    double lon = DEG2RAD(lon_deg);
    double lat = DEG2RAD(lat_deg);

    /* 卯酉圈曲率半径 N */
    double sin_lat = sin(lat);
    double cos_lat = cos(lat);
    double N = A / sqrt(1.0 - E2 * sin_lat * sin_lat);

    *X = (N + h) * cos_lat * cos(lon);
    *Y = (N + h) * cos_lat * sin(lon);
    *Z = (N * (1.0 - E2) + h) * sin_lat;
}

/* ECEF(XYZ) -> 大地坐标(经纬高)，采用迭代法 */
void ecef_to_geodetic_iterative(double X, double Y, double Z,
                                double *lon_deg, double *lat_deg, double *h) {
    double lon = atan2(Y, X);
    double p = sqrt(X * X + Y * Y);

    /* 纬度初值：常见近似 */
    double lat = atan2(Z, p * (1.0 - E2));
    double prev_lat;

    /* 迭代求解纬度 */
    const int max_iter = 100;
    const double tol = 1e-14; /* 弧度收敛阈值 */

    for (int i = 0; i < max_iter; ++i) {
        prev_lat = lat;
        double sin_lat = sin(lat);
        double N = A / sqrt(1.0 - E2 * sin_lat * sin_lat);
        lat = atan2(Z + E2 * N * sin_lat, p);

        if (fabs(lat - prev_lat) < tol) {
            break;
        }
    }

    double sin_lat = sin(lat);
    double N = A / sqrt(1.0 - E2 * sin_lat * sin_lat);
    *h = p / cos(lat) - N;

    *lon_deg = RAD2DEG(lon);
    *lat_deg = RAD2DEG(lat);
}

/* 在给定纬度，计算 1" 对应距离（纬度方向/经度方向） */
void one_arcsec_distance(double lat_deg, double *lat_meter, double *lon_meter) {
    double lat = DEG2RAD(lat_deg);
    double sin_lat = sin(lat);
    double cos_lat = cos(lat);

    /* 子午圈曲率半径 M 与卯酉圈曲率半径 N */
    double W = sqrt(1.0 - E2 * sin_lat * sin_lat);
    double N = A / W;
    double M = A * (1.0 - E2) / (W * W * W);

    /* 1" 对应弧度 */
    double arcsec = M_PI / 648000.0;

    *lat_meter = M * arcsec;          /* 纬度方向（子午线） */
    *lon_meter = N * cos_lat * arcsec;/* 经度方向（纬线） */
}

int main(void) {
    /* 已知坐标：经度 114°24'23.1455"，纬度 30°30'18.4323"，高程 20.258m */
    double lon_deg = dms_to_deg(114, 24, 23.1455);
    double lat_deg = dms_to_deg(30, 30, 18.4323);
    double h = 20.258;

    /* 1) 大地坐标 -> 空间直角坐标 */
    double X, Y, Z;
    geodetic_to_ecef(lon_deg, lat_deg, h, &X, &Y, &Z);

    /* 2) 空间直角坐标 -> 大地坐标（迭代） */
    double lon_back, lat_back, h_back;
    ecef_to_geodetic_iterative(X, Y, Z, &lon_back, &lat_back, &h_back);

    int lon_d, lon_m, lat_d, lat_m;
    double lon_s, lat_s;
    deg_to_dms(lon_back, &lon_d, &lon_m, &lon_s);
    deg_to_dms(lat_back, &lat_d, &lat_m, &lat_s);

    /* 3) 1" 对应距离 */
    double lat_1sec_m, lon_1sec_m;
    one_arcsec_distance(lat_deg, &lat_1sec_m, &lon_1sec_m);

    /* 4) 毫米精度下秒的小数位数估计 */
    double mm = 0.001;
    double sec_for_1mm_lat = mm / lat_1sec_m;
    double sec_for_1mm_lon = mm / lon_1sec_m;

    /* 向上估算需要的小数位数：10^-n <= sec_for_1mm */
    int n_lat = (int)ceil(-log10(sec_for_1mm_lat));
    int n_lon = (int)ceil(-log10(sec_for_1mm_lon));

    printf("已知大地坐标:\n");
    printf("  经度 = %.10f°\n", lon_deg);
    printf("  纬度 = %.10f°\n", lat_deg);
    printf("  高程 = %.3f m\n\n", h);

    printf("1) 空间直角坐标(ECEF):\n");
    printf("  X = %.4f m\n", X);
    printf("  Y = %.4f m\n", Y);
    printf("  Z = %.4f m\n\n", Z);

    printf("2) 由 XYZ 反算大地坐标(迭代):\n");
    printf("  经度 = %d°%d'%.8f\"  (%.10f°)\n", lon_d, lon_m, lon_s, lon_back);
    printf("  纬度 = %d°%d'%.8f\"  (%.10f°)\n", lat_d, lat_m, lat_s, lat_back);
    printf("  高程 = %.6f m\n\n", h_back);

    printf("3) 在该纬度 %.10f° 处，1\" 对应距离:\n", lat_deg);
    printf("  纬度方向(子午线): %.6f m\n", lat_1sec_m);
    printf("  经度方向(纬线)  : %.6f m\n\n", lon_1sec_m);

    printf("4) 若直角坐标精度约为 1 mm：\n");
    printf("  纬度秒值建议至少保留 %d 位小数 (1 mm ≈ %.8e\")\n", n_lat, sec_for_1mm_lat);
    printf("  经度秒值建议至少保留 %d 位小数 (1 mm ≈ %.8e\")\n", n_lon, sec_for_1mm_lon);
    printf("  综合建议：经纬度秒值保留 5 位小数较稳妥。\n");

    return 0;
}
