
from fractions import Fraction

# ===== các hàm hỗ trợ định dạng =====
def to_fraction(value):
    if isinstance(value, Fraction):
        return value
    return Fraction(str(value))

def dinh_dang_hien_thi(value):
    value = to_fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"

# In ma trận dạng đẹp đều
def in_ma_tran_dep(ma_tran):
    for dong in ma_tran:
        formatted_dong = [f"{dinh_dang_hien_thi(v):>3}" for v in dong]
        print("[ " + "  ".join(formatted_dong) + " ]")

# Tìm nghiệm tổng quát cho hệ có vô số nghiệm
def nghiem_tong_quat_he_vo_so_nghiem(ma_tran_sau_khu, b_sau_khu, pivot_rows, pivot_cols):
    so_dong_sau_khu = len(ma_tran_sau_khu)
    so_cot_sau_khu = len(ma_tran_sau_khu[0])
    # Tìm cột tự do (cột không phải cột pivot) 
    cot_tu_do_s = [col for col in range(so_cot_sau_khu) if col not in pivot_cols]
    # Đặt tên tham số: nếu 1 biến tự do thì là 't', nếu nhiều thì là 't1, t2, ...'
    tham_so = ["t" if len(cot_tu_do_s) == 1 else f"t{i + 1}" for i in range(len(cot_tu_do_s))]

    # Tính nghiệm riêng x0 
    x0 = [Fraction(0, 1)] * so_cot_sau_khu
    for idx in range(len(pivot_cols) - 1, -1, -1):
        rowPivot = pivot_rows[idx]  # Hàng chứa pivot
        colPivot = pivot_cols[idx]  # Cột chứa pivot
        tong_duoi_hang = Fraction(0, 1)
        # Tính tổng các phần tử bên phải pivot trong hàng này
        for j in range(colPivot + 1, so_cot_sau_khu):
            tong_duoi_hang += ma_tran_sau_khu[rowPivot][j] * x0[j]
        # Giải ngược: x[c] = (b - tổng_phía_phải) / a[r][c]
        x0[colPivot] = (b_sau_khu[rowPivot] - tong_duoi_hang) / ma_tran_sau_khu[rowPivot][colPivot]


    # Tính các vector cơ sở: mỗi vector ứng với một biến tự do
    co_so = []
    for cot_tu_do in cot_tu_do_s:
        # Vector v: cho biến tự do này = 1, các biến tự do khác = 0
        vec = [Fraction(0, 1)] * so_cot_sau_khu
        vec[cot_tu_do] = Fraction(1, 1)
        # Tính các thành phần còn lại bằng cách giải ngược từ dưới lên
        for idx in range(len(pivot_cols) - 1, -1, -1):
            rowPivot = pivot_rows[idx]
            colPivot = pivot_cols[idx]
            tong_duoi_hang = Fraction(0, 1)
            for j in range(colPivot + 1, so_cot_sau_khu):
                tong_duoi_hang += ma_tran_sau_khu[rowPivot][j] * vec[j]
            # Giải ngược: x[c] = -tổng_phía_phải / a[r][c]  
            vec[colPivot] = -tong_duoi_hang / ma_tran_sau_khu[rowPivot][colPivot]
        co_so.append(vec)


    # Xây dựng biểu thức: x_i = x0[i] + c1*v1[i]*t1 + c2*v2[i]*t2 + ...
    bieu_thuc_tung_an = []
    for i in range(so_cot_sau_khu):
        phan = []  
        hang_so = x0[i] 
        if hang_so != 0:
            phan.append(dinh_dang_hien_thi(hang_so))

        for j, vec in enumerate(co_so):
            he_so = vec[i]  
            if he_so == 0:
                continue  
            if he_so == 1:
                tmp = tham_so[j]
            elif he_so == -1:
                tmp = f"-{tham_so[j]}"
            else:
                tmp = f"{dinh_dang_hien_thi(he_so)}*{tham_so[j]}"
            phan.append(tmp)

        # Nối các phần lại, xử lý dấu âm
        s = " + ".join(phan).replace("+ -", "- ") if phan else "0"
        bieu_thuc_tung_an.append(f"x{i + 1} = {s}")

    # Trích xuất phần bên phải dấu '=' từ biểu thức để trả về
    nghiem_he = [dong.split("=", 1)[1].strip() for dong in bieu_thuc_tung_an]
    return nghiem_he

def print_gaussian_eliminate(A, b = None):
    # Lấy kết quả từ hàm gaussian_elimination để in ra ma trận đã khử, nghiệm và số lần đổi hàng
    ma_tran_sau_khu, nghiem_he, so_lan_doi_hang = gaussian_eliminate(A, b)
    he_vo_so_nghiem = isinstance(nghiem_he, list) and len(nghiem_he) > 0 and isinstance(nghiem_he[0], str)
    if nghiem_he is None:
        print("Hệ vô nghiệm.")
    elif he_vo_so_nghiem:
        print("Hệ có vô số nghiệm.")
    else:
        print("Hệ có nghiệm duy nhất.")

    print(f"Ma trận sau khi đã được khử:")
    # In ma trận dạng đẹp
    in_ma_tran_dep(ma_tran_sau_khu)
    # In kết quả tương ứng với loại nghiệm
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
# Thực hiện khử Gauss với phương pháp chọn pivot (partial pivoting)
def khu_ma_tran_ve_bac_thang(ma_tran_chuan_hoa, b_chuan_hoa):
    so_dong_chuan_hoa = len(ma_tran_chuan_hoa)      # Số hàng
    so_cot_chuan_hoa = len(ma_tran_chuan_hoa[0])   # Số cột

    so_lan_doi_hang = 0         # Đếm số lần hoán đổi hàng
    pivot_rows = [] # Lưu chỉ số hàng chứa pivot
    pivot_cols = [] # Lưu chỉ số cột chứa pivot
    k = 0           # Chỉ số hàng hiện tại (từ 0 trở đi)

    # Duyệt qua từng cột từ trái sang phải
    for col in range(so_cot_chuan_hoa):
        if k >= so_dong_chuan_hoa:  # Đã xử lý hết các hàng
            break

        # BỚC 1: TÌM PIVOT - tìm hàng có giá trị tuyệt đối lớn nhất ở cột này
        max_row = k
        for r in range(k + 1, so_dong_chuan_hoa):
            if abs(ma_tran_chuan_hoa[r][col]) > abs(ma_tran_chuan_hoa[max_row][col]):
                max_row = r

        # Nếu pivot = 0, cộ này không có pivot, chuyển sang cột tiếp theo
        if ma_tran_chuan_hoa[max_row][col] == 0:
            continue

        # BƯỚC 2: HOÁN HÀNG - đưa hàng có pivot lên vị trí k
        if max_row != k:
            # Hoán hàng của ma trận
            ma_tran_chuan_hoa[k], ma_tran_chuan_hoa[max_row] = ma_tran_chuan_hoa[max_row], ma_tran_chuan_hoa[k]
            # Hoán hàng của vector b tương ứng
            b_chuan_hoa[k], b_chuan_hoa[max_row] = b_chuan_hoa[max_row], b_chuan_hoa[k]
            so_lan_doi_hang += 1  # Tăng bộ đếm hoán hàng

        # BƯỚC 3: sắp xếp pivot
        pivot_rows.append(k)
        pivot_cols.append(col)

        # BƯỚC 4: khủ dưới - đặt các phần tử dưới pivot = 0
        for r in range(k + 1, so_dong_chuan_hoa):
            # Tính hệ số để khử: factor = a[r][col] / a[k][col]
            factor = ma_tran_chuan_hoa[r][col] / ma_tran_chuan_hoa[k][col]
            if factor == 0:
                # Nếu factor = 0 có nghĩa là a[r][col] = 0, không cần khử
                ma_tran_chuan_hoa[r][col] = Fraction(0, 1)
                continue

            # Thực hiện phép biến đổi hàng: hàng r = hàng r - factor * hàng k
            for c in range(col, so_cot_chuan_hoa):
                ma_tran_chuan_hoa[r][c] = ma_tran_chuan_hoa[r][c] - factor * ma_tran_chuan_hoa[k][c]
            b_chuan_hoa[r] = b_chuan_hoa[r] - factor * b_chuan_hoa[k]

        k += 1  # Chuyển sang hàng tiếp theo

    return ma_tran_chuan_hoa, b_chuan_hoa, so_lan_doi_hang, pivot_rows, pivot_cols


# Xác định loại nghiệm (vô nghiệm / vô số / duy nhất) và tính nghiệm nếu có
def tim_nghiem_he(ma_tran_sau_khu, b_sau_khu, pivot_rows, pivot_cols, A, b):
    so_dong_sau_khu = len(ma_tran_sau_khu)    # Số hàng
    so_cot_sau_khu = len(ma_tran_sau_khu[0]) # Số cột

    he_vo_nghiem = False # cờ flag để đánh dấu nếu hệ vô nghiệm
    for r in range(so_dong_sau_khu):
        # Kiểm tra tất cả phần tử trái = 0
        left_zero = all(ma_tran_sau_khu[r][c] == 0 for c in range(so_cot_sau_khu))
        # Nếu trái = 0 nhưng phải ≠ 0 thì vô nghiệm
        if left_zero and b_sau_khu[r] != 0:
            he_vo_nghiem = True
            break
    
    nghiem_he = []
    if he_vo_nghiem:
        # TRƯỜNG HỢP 1: HỆ VÔ NGHIỆM
        nghiem_he = None
    elif len(pivot_cols) < so_cot_sau_khu:
        # TRƯỜNG HỢP 2: HỆ CÓ VÔ SỐ NGHIỆM
        # (số pivot < số cột nghĩa là có biến tự do)
        nghiem_he = nghiem_tong_quat_he_vo_so_nghiem(ma_tran_sau_khu, b_sau_khu, pivot_rows, pivot_cols)
    else:
        # TRƯỜNG HỢP 3: HỆ CÓ NGHIỆM DUY NHẤT
        x = [Fraction(0, 1)] * so_cot_sau_khu
        # Duyệt từ pivot cuối cùng về trước
        for idx in range(len(pivot_cols) - 1, -1, -1):
            rowPivot = pivot_rows[idx]  # Hàng chứa pivot
            colPivot = pivot_cols[idx]  # Cột chứa pivot
            # Tính tổng các phần tử bên phải trong hàng
            tong_duoi_hang = Fraction(0, 1)
            for j in range(colPivot + 1, so_cot_sau_khu):
                tong_duoi_hang += ma_tran_sau_khu[rowPivot][j] * x[j]
            # Giải x[colPivot]: x[colPivot] = (b[rowPivot] - tổng) / a[rowPivot][colPivot]
            x[colPivot] = (b_sau_khu[rowPivot] - tong_duoi_hang) / ma_tran_sau_khu[rowPivot][colPivot]
        nghiem_he = x[:]

    return nghiem_he



# ===== HÀM CHÍNH - GIẢI HỆ PHƯƠNG TRÌNH TUYẾN TÍNH BẰNG GAUSSIAN ELIMINATION =====

# giải hệ Ax = b sử dụng phương pháp khử Gaussian
def gaussian_eliminate(A, b):
    # A và b phải là danh sách không rỗng
    if not A or not b:
        raise ValueError("A, b không được để trống.")
    
    # Kiểm tra A có phải là ma trận chữ nhật không
    if len(A) != len(b):
        raise ValueError("Số dòng của A phải bằng số dòng của b.")

    so_dong = len(A)  # Số phương trình
    so_cot = len(A[0])  # Số ẩn số

    # Kiểm tra tất cả hàng của A có cùng số cột không
    if any(len(row) != so_cot for row in A):
        raise ValueError("A phải là ma trận chữ nhật.")

    # CHUYỂN ĐỔI SANG FRACTION ĐỂ TÍNH TOÁN CHÍNH XÁC
    ma_tran_chuan_hoa = [[to_fraction(x) for x in row] for row in A]
    b_chuan_hoa = [to_fraction(x) for x in b]

    # BƯỚC 1: Khử Gaussian - đưa ma trận về dạng bậc thang
    ma_tran_sau_khu, b_sau_khu, so_lan_doi_hang, pivot_rows, pivot_cols = khu_ma_tran_ve_bac_thang(ma_tran_chuan_hoa, b_chuan_hoa)
    
    # BƯỚC 2: Tìm nghiệm dựa trên ma trận bậc thang
    nghiem_he = tim_nghiem_he(ma_tran_sau_khu, b_sau_khu, pivot_rows, pivot_cols, A, b)

    # trả về: (ma trận khử, nghiệm, số lần đổi hàng)
    return ma_tran_sau_khu, nghiem_he, so_lan_doi_hang


# ===== HÀM IN KẾT QUẢ BACK SUBSTITUTION =====

def print_back_substitution(U, c):
    # lưu kết quả back substitution để in ra sau khi đã phân loại nghiệm
    nghiem_he = back_substitution(U, c)
    print(f"Ma trận U: ")
    # In ma trận dạng đẹp
    in_ma_tran_dep(U)
    print(f"Vector c: {[dinh_dang_hien_thi(v) for v in c]}")
    
    if nghiem_he is None:
        print("Nghiệm hệ: vô nghiệm")
    elif isinstance(nghiem_he, list) and len(nghiem_he) > 0 and isinstance(nghiem_he[0], str):
        print(f"Nghiệm hệ: {nghiem_he}")
    elif isinstance(nghiem_he, list):
        print(f"Nghiệm hệ: {[dinh_dang_hien_thi(v) for v in nghiem_he]}")
    print()
# ===== GIẢI NGƯỢC - BACK SUBSTITUTION CHO MA TRẬN TAM GIÁC TRÊN =====
def back_substitution(U, c):
    if not U or not c:
        raise ValueError("U, c không được để trống.")

    so_dong = len(U) # Số hàng của ma trận tam giác trên
    so_cot = len(U[0]) # Số cột của ma trận tam giác trên
    if len(c) != so_dong:
        raise ValueError("Số phần tử của c phải bằng số hàng của U.")
    if so_dong != so_cot:
        raise ValueError("U phải là ma trận vuông (n x n).")

    U_chuan_hoa = [[to_fraction(val) for val in row] for row in U] 
    c_chuan_hoa = [to_fraction(val) for val in c]
    
    # Kiểm tra nếu có hàng 0 = 0 0 ... 0 | c (c ≠ 0) thì hệ vô nghiệm
    for i in range(so_dong):
        left_zero = all(U_chuan_hoa[i][j] == 0 for j in range(so_cot))
        if left_zero and c_chuan_hoa[i] != 0:
            return None  # Vô nghiệm
    
    pivot_rows = []
    pivot_cols = []
    
    # Xác định vị trí pivot
    for i in range(so_dong):
        if U_chuan_hoa[i][i] != 0:
            pivot_rows.append(i)
            pivot_cols.append(i)
    
    # sử dụng hàm tim_nghiem_he để phân loại và tìm nghiệm (nếu có) cho hệ Ux = c
    nghiem_he = tim_nghiem_he(U_chuan_hoa, c_chuan_hoa, pivot_rows, pivot_cols, U, c)
    
    return nghiem_he


# ===== HÀM MAIN =====

def main():
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


if __name__ == "__main__":
    main()