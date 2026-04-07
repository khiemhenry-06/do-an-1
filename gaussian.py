import random
import time

EPS = 1e-12


def to_float(value):
    # Chuẩn hóa dữ liệu đầu vào về float.
    return float(value)

# ===== các hàm hỗ trợ định dạng =====
def dinh_dang_hien_thi(value, cham_thap_phan=6):
    if isinstance(value, str):
        return value
    if abs(value) < EPS:
        value = 0.0
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.{cham_thap_phan}g}"


def in_ma_tran_dep(ma_tran):
    # In ma trận dạng đẹp, căn lề để kết quả dễ quan sát.
    for dong in ma_tran:
        formatted_dong = [f"{dinh_dang_hien_thi(v):>10}" for v in dong]
        print("[ " + "  ".join(formatted_dong) + " ]")


# Tìm nghiệm tổng quát cho hệ có vô số nghiệm
def nghiem_tong_quat_he_vo_so_nghiem(ma_tran_sau_khu, b_sau_khu, pivot_rows, pivot_cols, eps=EPS):
    so_cot_sau_khu = len(ma_tran_sau_khu[0])
    # Tìm các cột tự do (không phải cột pivot).
    cot_tu_do_s = [col for col in range(so_cot_sau_khu) if col not in pivot_cols]
    # Đặt tên tham số: 1 biến tự do -> t, nhiều biến -> t1, t2, ...
    tham_so = ["t" if len(cot_tu_do_s) == 1 else f"t{i + 1}" for i in range(len(cot_tu_do_s))]

    x0 = [0.0] * so_cot_sau_khu

    # Tính nghiệm đặc biệt x0 với các tham số tự do bằng 0
    for idx in range(len(pivot_cols) - 1, -1, -1):
        row_pivot = pivot_rows[idx]
        col_pivot = pivot_cols[idx]
        tong_duoi_hang = 0.0
        for j in range(col_pivot + 1, so_cot_sau_khu):
            tong_duoi_hang += ma_tran_sau_khu[row_pivot][j] * x0[j]
        x0[col_pivot] = (b_sau_khu[row_pivot] - tong_duoi_hang) / ma_tran_sau_khu[row_pivot][col_pivot]

    co_so = []
    for cot_tu_do in cot_tu_do_s:
        # Tạo vector cơ sở: biến tự do hiện tại = 1, các biến tự do khác = 0.
        vec = [0.0] * so_cot_sau_khu
        vec[cot_tu_do] = 1.0
        for idx in range(len(pivot_cols) - 1, -1, -1):
            row_pivot = pivot_rows[idx]
            col_pivot = pivot_cols[idx]
            tong_duoi_hang = 0.0
            for j in range(col_pivot + 1, so_cot_sau_khu):
                tong_duoi_hang += ma_tran_sau_khu[row_pivot][j] * vec[j]
            vec[col_pivot] = -tong_duoi_hang / ma_tran_sau_khu[row_pivot][col_pivot]
        co_so.append(vec)

    bieu_thuc_tung_an = []
    for i in range(so_cot_sau_khu):
        # Xây dựng x_i = hằng số + tổ hợp tuyến tính của các tham số.
        phan = []
        hang_so = x0[i]
        if abs(hang_so) > eps:
            phan.append(dinh_dang_hien_thi(hang_so))

        for j, vec in enumerate(co_so):
            he_so = vec[i]
            if abs(he_so) <= eps:
                continue
            if abs(he_so - 1.0) <= eps:
                tmp = tham_so[j]
            elif abs(he_so + 1.0) <= eps:
                tmp = f"-{tham_so[j]}"
            else:
                tmp = f"{dinh_dang_hien_thi(he_so)}*{tham_so[j]}"
            phan.append(tmp)

        s = " + ".join(phan).replace("+ -", "- ") if phan else "0"
        bieu_thuc_tung_an.append(f"x{i + 1} = {s}")

    nghiem_he = [dong.split("=", 1)[1].strip() for dong in bieu_thuc_tung_an]
    return nghiem_he


def print_gaussian_eliminate(A, b=None, eps=EPS):
    # Lấy kết quả từ hàm gaussian_eliminate để in đồng bộ.
    ma_tran_sau_khu, nghiem_he, so_lan_doi_hang = gaussian_eliminate(A, b, eps=eps)
    he_vo_so_nghiem = isinstance(nghiem_he, list) and len(nghiem_he) > 0 and isinstance(nghiem_he[0], str)
    if nghiem_he is None:
        print("Hệ vô nghiệm.")
    elif he_vo_so_nghiem:
        print("Hệ có vô số nghiệm.")
    else:
        print("Hệ có nghiệm duy nhất.")

    print("Ma trận sau khi đã được khử:")
    in_ma_tran_dep(ma_tran_sau_khu)

    if nghiem_he is None:
        print("Nghiệm hệ: vô nghiệm")
    elif he_vo_so_nghiem:
        print(f"Nghiệm hệ: {nghiem_he}")
    elif nghiem_he:
        print(f"Nghiệm hệ: {[dinh_dang_hien_thi(v) for v in nghiem_he]}")
    else:
        print("Nghiệm hệ: []")
    print(f"Số lần hoán đổi hàng: {so_lan_doi_hang}")
    print()


# Thực hiện khử Gauss với partial pivoting
def khu_ma_tran_ve_bac_thang(ma_tran_chuan_hoa, b_chuan_hoa, eps=EPS):
    # so_dong_chuan_hoa: số phương trình, so_cot_chuan_hoa: số ẩn.
    so_dong_chuan_hoa = len(ma_tran_chuan_hoa)
    so_cot_chuan_hoa = len(ma_tran_chuan_hoa[0])

    so_lan_doi_hang = 0
    pivot_rows = []
    pivot_cols = []
    k = 0

    for col in range(so_cot_chuan_hoa):
        if k >= so_dong_chuan_hoa:
            break

        # BƯỚC 1: tìm pivot có trị tuyệt đối lớn nhất trên cột hiện tại.
        max_row = k
        for r in range(k + 1, so_dong_chuan_hoa):
            if abs(ma_tran_chuan_hoa[r][col]) > abs(ma_tran_chuan_hoa[max_row][col]):
                max_row = r

        # Nếu cột này không có pivot hợp lệ thì bỏ qua.
        if abs(ma_tran_chuan_hoa[max_row][col]) <= eps:
            continue

        # BƯỚC 2: hoán đổi hàng để đưa pivot về vị trí [k][col].
        if max_row != k:
            ma_tran_chuan_hoa[k], ma_tran_chuan_hoa[max_row] = ma_tran_chuan_hoa[max_row], ma_tran_chuan_hoa[k]
            b_chuan_hoa[k], b_chuan_hoa[max_row] = b_chuan_hoa[max_row], b_chuan_hoa[k]
            so_lan_doi_hang += 1

        # Lưu vị trí pivot.
        pivot_rows.append(k)
        pivot_cols.append(col)

        # BƯỚC 3: khử các phần tử bên dưới pivot.
        for r in range(k + 1, so_dong_chuan_hoa):
            factor = ma_tran_chuan_hoa[r][col] / ma_tran_chuan_hoa[k][col]
            if abs(factor) <= eps:
                ma_tran_chuan_hoa[r][col] = 0.0
                continue

            for c in range(col, so_cot_chuan_hoa):
                ma_tran_chuan_hoa[r][c] = ma_tran_chuan_hoa[r][c] - factor * ma_tran_chuan_hoa[k][c]
                if abs(ma_tran_chuan_hoa[r][c]) <= eps:
                    ma_tran_chuan_hoa[r][c] = 0.0

            b_chuan_hoa[r] = b_chuan_hoa[r] - factor * b_chuan_hoa[k]
            if abs(b_chuan_hoa[r]) <= eps:
                b_chuan_hoa[r] = 0.0

        k += 1

    # Trả về ma trận bậc thang, vector b đã biến đổi và thông tin pivot.
    return ma_tran_chuan_hoa, b_chuan_hoa, so_lan_doi_hang, pivot_rows, pivot_cols


def tim_nghiem_he(ma_tran_sau_khu, b_sau_khu, pivot_rows, pivot_cols, A, b, eps=EPS):
    so_dong_sau_khu = len(ma_tran_sau_khu)
    so_cot_sau_khu = len(ma_tran_sau_khu[0])

    he_vo_nghiem = False
    for r in range(so_dong_sau_khu):
        # Nếu vế trái là hàng 0, vế phải khác 0 -> hệ vô nghiệm.
        left_zero = all(abs(ma_tran_sau_khu[r][c]) <= eps for c in range(so_cot_sau_khu))
        if left_zero and abs(b_sau_khu[r]) > eps:
            he_vo_nghiem = True
            break

    # Phân loại: vô nghiệm / vô số nghiệm / nghiệm duy nhất.
    if he_vo_nghiem:
        nghiem_he = None
    elif len(pivot_cols) < so_cot_sau_khu:
        nghiem_he = nghiem_tong_quat_he_vo_so_nghiem(ma_tran_sau_khu, b_sau_khu, pivot_rows, pivot_cols, eps=eps)
    else:
        # TRƯỜNG HỢP nghiệm duy nhất: thế ngược từ dưới lên.
        x = [0.0] * so_cot_sau_khu
        for idx in range(len(pivot_cols) - 1, -1, -1):
            row_pivot = pivot_rows[idx]
            col_pivot = pivot_cols[idx]
            tong_duoi_hang = 0.0
            for j in range(col_pivot + 1, so_cot_sau_khu):
                tong_duoi_hang += ma_tran_sau_khu[row_pivot][j] * x[j]
            x[col_pivot] = (b_sau_khu[row_pivot] - tong_duoi_hang) / ma_tran_sau_khu[row_pivot][col_pivot]
        nghiem_he = x[:]

    return nghiem_he


def gaussian_eliminate(A, b, eps=EPS):
    # A, b phải không rỗng.
    if not A or not b:
        raise ValueError("A, b không được để trống.")

    # Số dòng của A phải bằng số phần tử của b.
    if len(A) != len(b):
        raise ValueError("Số dòng của A phải bằng số dòng của b.")

    # Kiểm tra A là ma trận chữ nhật.
    so_cot = len(A[0])
    if any(len(row) != so_cot for row in A):
        raise ValueError("A phải là ma trận chữ nhật.")

    # Chuẩn hóa sang float để tính nhanh.
    ma_tran_chuan_hoa = [[to_float(x) for x in row] for row in A]
    b_chuan_hoa = [to_float(x) for x in b]

    # BƯỚC 1: Khử Gaussian đưa về bậc thang.
    ma_tran_sau_khu, b_sau_khu, so_lan_doi_hang, pivot_rows, pivot_cols = khu_ma_tran_ve_bac_thang(
        ma_tran_chuan_hoa, b_chuan_hoa, eps=eps
    )

    # BƯỚC 2: Tìm nghiệm từ ma trận bậc thang.
    nghiem_he = tim_nghiem_he(ma_tran_sau_khu, b_sau_khu, pivot_rows, pivot_cols, A, b, eps=eps)
    return ma_tran_sau_khu, nghiem_he, so_lan_doi_hang


def print_back_substitution(U, c, eps=EPS):
    # In kết quả back substitution theo định dạng thống nhất.
    nghiem_he = back_substitution(U, c, eps=eps)
    print("Ma trận U:")
    in_ma_tran_dep(U)
    print(f"Vector c: {[dinh_dang_hien_thi(v) for v in c]}")

    if nghiem_he is None:
        print("Nghiệm hệ: vô nghiệm")
    elif isinstance(nghiem_he, list) and len(nghiem_he) > 0 and isinstance(nghiem_he[0], str):
        print(f"Nghiệm hệ: {nghiem_he}")
    elif isinstance(nghiem_he, list):
        print(f"Nghiệm hệ: {[dinh_dang_hien_thi(v) for v in nghiem_he]}")
    print()


def back_substitution(U, c, eps=EPS):
    # U, c phải không rỗng.
    if not U or not c:
        raise ValueError("U, c không được để trống.")

    # U phải là ma trận vuông và kích thước hợp với c.
    so_dong = len(U)
    so_cot = len(U[0])
    if len(c) != so_dong:
        raise ValueError("Số phần tử của c phải bằng số hàng của U.")
    if so_dong != so_cot:
        raise ValueError("U phải là ma trận vuông (n x n).")

    U_chuan_hoa = [[to_float(val) for val in row] for row in U]
    c_chuan_hoa = [to_float(val) for val in c]

    # Kiểm tra hàng 0 | c != 0 -> vô nghiệm.
    for i in range(so_dong):
        left_zero = all(abs(U_chuan_hoa[i][j]) <= eps for j in range(so_cot))
        if left_zero and abs(c_chuan_hoa[i]) > eps:
            return None

    # Xác định các pivot trên đường chéo của U.
    pivot_rows = []
    pivot_cols = []
    for i in range(so_dong):
        if abs(U_chuan_hoa[i][i]) > eps:
            pivot_rows.append(i)
            pivot_cols.append(i)

    # Tái sử dụng hàm phân loại nghiệm chung cho Ux = c.
    nghiem_he = tim_nghiem_he(U_chuan_hoa, c_chuan_hoa, pivot_rows, pivot_cols, U, c, eps=eps)
    return nghiem_he


def random_matrix_vector(n, m):
    A = [[random.randint(-10, 10) for _ in range(m)] for _ in range(n)]
    b = [random.randint(-10, 10) for _ in range(n)]
    return A, b


def main_test():
    print("\n" + "="*50)
    print("TEST HÀM GAUSSIAN ELIMINATION")
    print("="*50 + "\n")
    # TEST 1: Hệ có nghiệm duy nhất 
    A1 = [[1, 2, 0, 2], [3, 5, -1, 6], [2, 4, 1, 2], [2, 0, -7, 11]]
    b1 = [6, 17, 12, 7]
    print_gaussian_eliminate(A1, b1)

    # TEST 2: Hệ vô nghiệm 
    A2 = [[2, -4, -1], [1, -3, 1], [3, -5, -3]]
    b2 = [1, 1, 2]
    print_gaussian_eliminate(A2, b2)

    # TEST 3: Hệ có vô số nghiệm 
    A3 = [[1, -2, -1], [2, -3, 1], [3, -5, 0], [1, 0, 5]]
    b3 = [1, 6, 7, 9]
    print_gaussian_eliminate(A3, b3)

    # ===== TEST HÀM BACK_SUBSTITUTION =====
    print("\n" + "="*50)
    print("TEST HÀM BACK_SUBSTITUTION")
    print("="*50 + "\n")
    
    # TEST 1: Ma trận tam giác trên - Hệ có nghiệm duy nhất
    U1 = [[2, 1, 3], [0, 4, 2], [0, 0, 1]]
    c1 = [13, 18, 5]
    print_back_substitution(U1, c1)
    
    # TEST 2: Ma trận tam giác trên - Hệ vô nghiệm
    U2 = [[1, 2], [0, 0]]
    c2 = [5, 1]
    print_back_substitution(U2, c2)
    
    # TEST 3: Ma trận tam giác trên - Hệ vô số nghiệm
    U3 = [[1, 1, 2], [0, 0, 0], [0, 0, 0]]
    c3 = [5, 0, 0]
    print_back_substitution(U3, c3)


def main_test_time():
    sizes = [10, 20, 50, 100, 200, 500, 1000]
    for size in sizes:
        A, b = random_matrix_vector(size, size)
        start_time = time.time()
        gaussian_eliminate(A, b)
        end_time = time.time()
        print(f"Thời gian giải hệ {size}x{size}: {end_time - start_time:.4f} giây")


if __name__ == "__main__":
    main_test()
    #main_test_time()
