from fractions import Fraction
from gaussian import *

# ===== TÍNH MA TRẬN NGHỊCH ĐẢO BẰNG PHƯƠNG PHÁP GAUSS-JORDAN =====
def inverse(A):
    # Kiểm tra A có rỗng ko
    if not A:
        raise ValueError("A không được để trống.")
    
    # Kiểm tra A có vuông ko
    so_dong = len(A)
    so_cot = len(A[0])
    if so_dong != so_cot:
        raise ValueError("A phải là ma trận vuông.")

    ma_tran_chuan_hoa = [[to_fraction(x) for x in row] for row in A]
    
    # tạo ma trận đơn vị I
    ma_tran_I = [[Fraction(1, 1) if i == j else Fraction(0, 1) for j in range(so_cot)] for i in range(so_dong)]
    
    # tạo ma trận [A | I]
    ma_tran_tang_cuong = [ma_tran_chuan_hoa[i] + ma_tran_I[i] for i in range(so_dong)]

    # GAUSS-JORDAN
    # BƯỚC 1: Khử xuống 
    for col in range(so_cot):
        # Tìm pivot - hàng có giá trị tuyệt đối lớn nhất
        max_row = col
        for row in range(col + 1, so_dong):
            if abs(ma_tran_tang_cuong[row][col]) > abs(ma_tran_tang_cuong[max_row][col]):
                max_row = row
        
        # Nếu pivot = 0, ma trận không khả nghịch
        if ma_tran_tang_cuong[max_row][col] == 0:
            return None
        
        # đổi hàng nếu cần
        if max_row != col:
            ma_tran_tang_cuong[col], ma_tran_tang_cuong[max_row] = ma_tran_tang_cuong[max_row], ma_tran_tang_cuong[col]
        
        # Khử các phần tử dưới pivot
        pivot = ma_tran_tang_cuong[col][col]
        for row in range(col + 1, so_dong):
            if ma_tran_tang_cuong[row][col] != 0:
                factor = ma_tran_tang_cuong[row][col] / pivot
                for c in range(2 * so_cot):
                    ma_tran_tang_cuong[row][c] = ma_tran_tang_cuong[row][c] - factor * ma_tran_tang_cuong[col][c]
    
    # BƯỚC 2: Khử lên 
    for col in range(so_cot - 1, -1, -1):
        # Chuẩn hóa hàng col để pivot = 1
        pivot = ma_tran_tang_cuong[col][col]
        if pivot == 0:
            return None  # Ma trận không khả nghịch
        
        for c in range(2 * so_cot):
            ma_tran_tang_cuong[col][c] = ma_tran_tang_cuong[col][c] / pivot
        
        # Khử các phần tử phía trên pivot
        for row in range(col):
            if ma_tran_tang_cuong[row][col] != 0:
                factor = ma_tran_tang_cuong[row][col]
                for c in range(2 * so_cot):
                    ma_tran_tang_cuong[row][c] = ma_tran_tang_cuong[row][c] - factor * ma_tran_tang_cuong[col][c]

    # Lấy má trận nghịch đảo (phần bên phải)
    A_nghich_dao = [row[so_cot:2 * so_cot] for row in ma_tran_tang_cuong]
    
    return A_nghich_dao


def print_inverse(A):
    """In ma trận A và ma trận nghịch đảo A^-1"""
    print(f"Ma trận A:")
    in_ma_tran_dep(A)
    print()
    
    A_nghich_dao = inverse(A)
    
    if A_nghich_dao is None:
        print("Ma trận A không khả nghịch do det(A) = 0.\n")
    else:
        print(f"Ma trận nghịch đảo A^-1:")
        in_ma_tran_dep(A_nghich_dao)
        print()


# ===== HÀM MAIN =====

def main():
    print("\n" + "="*50)
    print("TEST HÀM INVERSE (MA TRẬN NGHỊCH ĐẢO)")
    print("="*50 + "\n")
    
    # TEST 1: Ma trận 2x2
    A1 = [[2, 3], [1, 5]]
    print("TEST 1: Ma trận 2x2")
    print_inverse(A1)
    
    # TEST 2: Ma trận 3x3
    A2 = [[1, 2, 3], [0, 1, 4], [5, 6, 0]]
    print("TEST 2: Ma trận 3x3")
    print_inverse(A2)
    
    # TEST 3: Ma trận đơn vị (Identity) - nghịch đảo của I chính là I
    A3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    print("TEST 3: Ma trận đơn vị 3x3")
    print_inverse(A3)
    
    # TEST 4: Ma trận không khả nghịch (singular)
    A4 = [[1, 2, 3], [2, 4, 6], [0, 1, 2]]
    print("TEST 4: Ma trận không khả nghịch")
    print_inverse(A4)
    
    # TEST 5: Ma trận 4x4
    A5 = [[1, 0, 2, -1], [3, 0, 0, 5], [2, 1, 4, -3], [1, 0, 5, 0]]
    print("TEST 5: Ma trận 4x4")
    print_inverse(A5)


if __name__ == "__main__":
    main()

    