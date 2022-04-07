# Hàm ghép trơn (SPLINE) - Nhóm 22 

## Input, ouput
* Input: Bộ điểm xi,yi sẽ nhập vào file input.txt
* Output: các đa thức Sline bậc 3,bậc 1,bậc 2  + vẽ đa các thức, chấm các điểm trong bộ điểm

## Hướng dẫn sử dụng
Chạy như file python bình thường.

## Phân tích ưu nhược điểm thuật toán
* Ưu điểm: 
    * Sử dụng trong bộ dữ liệu lớn, nhiều mốc nội suy. 
        Khi tăng thêm các mốc, ta vẫn thu được các đa thức với bậc cố định 
* Nhược điểm: 
    * Phụ thuộc vào việc chọn mốc nội suy
        Nếu ta chọn mốc nội suy hợp lý thì đa thức thu được sẽ chính xác
        Ngược lại, đa thức thu được có sai số khá lớn