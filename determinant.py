from gaussian import *

# ===== TÍNH ĐỊNH THỨC BẰNG KHỬ GAUSS =====
def determinant(A):
    # Kiểm tra A rỗng hay ko
    if not A:
        raise ValueError("A không được để trống.")
    
    # Kiểm tra ma trận vuông
    so_dong = len(A)
    so_cot = len(A[0])
    if so_dong != so_cot:
        raise ValueError("A phải là ma trận vuông.")

    ma_tran_chuan_hoa = [[to_float(x) for x in row] for row in A]
    b_chuan_hoa = [0.0] * so_dong

    ma_tran_sau_khu, b_sau_khu, so_lan_doi_hang, pivot_rows, pivot_cols = khu_ma_tran_ve_bac_thang(ma_tran_chuan_hoa, b_chuan_hoa, eps=EPS)

    # Kiểm tra nếu có hàng bằng 0 thì det = 0
    for i in range(so_dong):
        if all(abs(ma_tran_sau_khu[i][j]) <= EPS for j in range(so_cot)):
            return 0.0

    # Tính tích các phần tử trên đường chéo
    det_product = 1.0
    for i in range(so_cot):
        det_product *= ma_tran_sau_khu[i][i]

    # Định thức = (-1)^(số lần hoán hàng) × tích các phần tử trên đường chéo
    det = det_product * ((-1.0) ** so_lan_doi_hang)

    return det

def print_determinant(A):
    det = determinant(A)
    print(f"Ma trận: ")
    # In ma trận dạng đẹp
    in_ma_tran_dep(A)
    print(f"Định thức: {dinh_dang_hien_thi(det)}\n")    
# ===== HÀM MAIN =====

def main():
    print("\n" + "="*50)
    print("TEST HÀM DETERMINANT (ĐỊNH THỨC)")
    print("="*50 + "\n")
    
    # TEST 1: Ma trận 2x2
    A1 = [[2, 3], [1, 5]]
    print_determinant(A1)
    
    # TEST 2: Ma trận 3x3 
    A2 = [[1, 2, 3], [0, 1, 4], [5, 6, 0]]
    print_determinant(A2)
    
    # TEST 3: Ma trận có các hàng phụ thuộc tuyến tính (hàng 2 = 2*hàng 1) - Định thức = 0
    A3 = [[1, 2, 3], [2, 4, 6], [0, 1, 2]]
    print_determinant(A3)
    
    # TEST 4: Ma trận 4x4
    A4 = [[1, 0, 2, -1], [3, 0, 0, 5], [2, 1, 4, -3], [1, 0, 5, 0]]
    print_determinant(A4)
    
    # TEST 5: Ma trận đơn vị - Định thức = 1
    A5 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    print_determinant(A5)
    
    # TEST 6: Ma trận đường chéo - Định thức = tích các phần tử đường chéo
    A6 = [[2, 0, 0, 0], [0, 3, 0, 0], [0, 0, -1, 0], [0, 0, 0, 4]]
    print_determinant(A6)


if __name__ == "__main__":
    main()