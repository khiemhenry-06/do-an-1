# Hướng dẫn lấy giá trị trả về để dùng tiếp

Tài liệu này giải thích nhanh cách dùng các hàm trong `gaussian.py`, đặc biệt là hàm `gaussian_eliminate(A, b)`.

## 1) Hàm `gaussian_eliminate(A, b)` trả về gì?

Hàm trả về một tuple gồm 3 phần tử:

- phần tử `[0]`: `ma_tran_sau_khu` (ma trận dạng bậc thang, kiểu `float`)
- phần tử `[1]`: `nghiem_he`
- phần tử `[2]`: `so_lan_hoan_doi_hang` (`int`)

Dạng tổng quát:

```python
ma_tran_sau_khu, nghiem_he, so_lan_doi_hang = gaussian_eliminate(A, b)
```

## 2) Ý nghĩa của `nghiem_he` (phần tử `[1]`)

`nghiem_he` có 3 trường hợp:

- `None`: hệ vô nghiệm
- `list[float]`: hệ có nghiệm duy nhất
- `list[str]`: hệ có vô số nghiệm (dạng tham số)

Ví dụ kiểm tra nhanh:

```python
if nghiem_he is None:
    print("Vô nghiệm")
elif len(nghiem_he) > 0 and isinstance(nghiem_he[0], str):
    print("Vô số nghiệm, dạng tham số:", nghiem_he)
else:
    print("Nghiệm duy nhất:", nghiem_he)
```

## 3) Cách xử lý sai số số thực

Dùng `float`, nên khi so sánh với 0 sẽ dùng ngưỡng `EPS`:

```python
EPS = 1e-12

if abs(x) <= EPS:
    # coi như bằng 0
    ...
```

Bạn có thể truyền `int`/`float` bình thường. Hàm sẽ tự chuyển về `float` bên trong.

## 4) Khác nhau giữa `gaussian_eliminate` và `print_gaussian_eliminate`?

- `gaussian_eliminate(A, b)`: dùng để lấy dữ liệu trả về và xử lý tiếp trong code.
- `print_gaussian_eliminate(A, b)`: dùng để in kết quả ra màn hình.

Lưu ý: hiện tại `print_gaussian_eliminate` chỉ in, không `return` nghiệm để dùng tiếp. Nếu cần xử lý tiếp, hãy gọi `gaussian_eliminate`.

## 5) Mẫu dùng để “lấy rồi xài” (để gọi cho phần khác)

```python
A = [[2, 3], [1, 5]]
b = [7, 11]

R, x, swap_count = gaussian_eliminate(A, b)

if x is None:
    # xử lý vô nghiệm
    pass
elif len(x) > 0 and isinstance(x[0], str):
    # xử lý vô số nghiệm (tham số)
    pass
else:
    # xử lý nghiệm duy nhất
    # x là list float, dùng tiếp được ngay
    pass
```
