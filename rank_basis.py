from fractions import Fraction
from gaussian import *


def rank_and_basis(A):
	if not A:
		raise ValueError("A không được để trống.")

	so_dong = len(A)
	so_cot = len(A[0])
	if any(len(row) != so_cot for row in A):
		raise ValueError("A phải là ma trận chữ nhật.")

	A_chuan_hoa = [[to_fraction(x) for x in row] for row in A]
	b_chuan_hoa = [Fraction(0, 1)] * so_dong

	# Sao chép ma trận A để khử mà không làm thay đổi A gốc, vì cần A gốc để lấy cơ sở không gian cột
	ma_tran_sao = [row[:] for row in A_chuan_hoa]    
	ket_qua_khu = khu_ma_tran_ve_bac_thang(ma_tran_sao, b_chuan_hoa)
	R = ket_qua_khu[0]          # Ma trận bậc thang
	pivot_rows = ket_qua_khu[3] # Các hàng pivot
	pivot_cols = ket_qua_khu[4] # Các cột pivot

	rank = len(pivot_cols)      # Hạng bằng số cột pivot

	# Cơ sở không gian cột: lấy các cột pivot từ ma trận gốc
	khong_gian_cot = []
	for col in pivot_cols:
		khong_gian_cot.append([A_chuan_hoa[row][col] for row in range(so_dong)])

	# Cơ sở không gian dòng: các hàng khác 0 của ma trận bậc thang
	khong_gian_dong = []
	for row in range(so_dong):
		if any(R[row][col] != 0 for col in range(so_cot)):
			khong_gian_dong.append(R[row][:])

	# Cơ sở không gian nghiệm: dựng theo biến tự do
	cot_tu_do_s = [colPivot for colPivot in range(so_cot) if colPivot not in pivot_cols]
	khong_gian_nghiem= []
	for cot_tu_do in cot_tu_do_s:
		vec = [Fraction(0, 1)] * so_cot	
		vec[cot_tu_do] = Fraction(1, 1)

		# Giải ngược cho các biến pivot
		for idx in range(len(pivot_cols) - 1, -1, -1):
			rowPivot = pivot_rows[idx]
			colPivot = pivot_cols[idx]
			tong = Fraction(0, 1)
			for j in range(colPivot + 1, so_cot):
				tong += R[rowPivot][j] * vec[j]
			vec[colPivot] = -tong / R[rowPivot][colPivot]

		khong_gian_nghiem.append(vec)

	return rank, khong_gian_cot, khong_gian_dong, khong_gian_nghiem



def print_rank_and_basis(A):
	rank, khong_gian_cot, khong_gian_dong, khong_gian_nghiem = rank_and_basis(A)
	print("Ma trận:")
	in_ma_tran_dep(A)
	print(f"Hạng (rank): {rank}\n")

	print("Cơ sở không gian cột:")
	if not khong_gian_cot:
		print("  []")
	else:
		for vec in khong_gian_cot:
			print("  [" + " ".join(dinh_dang_hien_thi(v) for v in vec) + "]")
	print()

	print("Cơ sở không gian dòng:")
	if not khong_gian_dong:
		print("  []")
	else:
		for vec in khong_gian_dong:
			print("  [" + " ".join(dinh_dang_hien_thi(v) for v in vec) + "]")
	print()

	print("Cơ sở không gian nghiệm:")
	if not khong_gian_nghiem:
		print("  Không gian nghiệm chỉ chứa vector 0.")
	else:
		for vec in khong_gian_nghiem:
			print("  [" + " ".join(dinh_dang_hien_thi(v) for v in vec) + "]")
	print()

def main():
	print("\n" + "=" * 50)
	print("TEST HÀM RANK_AND_BASIS")
	print("=" * 50 + "\n")

	# Test 1: rank đầy đủ, null space rỗng
	A1 = [[1, 2], [3, 4]]
	print("TEST 1: Ma trận vuông khả nghịch")
	print_rank_and_basis(A1)

	# Test 2: rank thiếu, có null space
	A2 = [[1, 2, 3], [2, 4, 6], [0, 1, 2]]
	print("TEST 2: Ma trận phụ thuộc tuyến tính")
	print_rank_and_basis(A2)

	# Test 3: ma trận chữ nhật m < n
	A3 = [[1, 0, 2, -1], [2, 1, 4, -3]]
	print("TEST 3: Ma trận chữ nhật")
	print_rank_and_basis(A3)


if __name__ == "__main__":
	main()

