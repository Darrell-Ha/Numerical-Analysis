/**
 *  ?           <<<< Nội suy ngược >>>>
 *  * Dùng khi cần giải f(x) = y nhưng f(x) cho dưới dạng các giá trị rời rạc
 *  ! Nguyên tắc:
 *      * 1. Giải thô: Xác định khoảng cách ly nghiệm
 *          ? Khoảng đơn điệu ?
 *          ? Khoảng có chứa nghiệm ?
 *      * 2. Tìm nghiệm trên khoảng cách ly nghiệm tìm được
 *              ! Cố gắng chọn khoảng cách ly không nhỏ 
 *      *    Hoặc khoảng không cách ly nghiệm nhưng có nghiệm duy nhất
 *  ##################################################################### 
 *  ? Tổng quát:
 *  TODO-1: Tìm các khoảng là song ánh => return danh sách điểm đầu - cuối (left, right) của khoảng song ánh và chứa nghiệm
 *  TODO-2: Duyệt các khoảng song ánh tìm được => return danh sách (left, right) của các khoảng song ánh tương đương để tiến hành nội suy trên đó
 *  TODO-3: Duyệt các khoảng thu được ở TODO-2 => Mỗi khoảng cho một nghiệm => danh sách các nghiệm xấp xỉ
 * 
 *  ? Cụ thể:
 *  TODO-1: Tìm các khoảng là song ánh (tìm song ánh để đảm bảo điều kiện bài toàn nội suy đa thức)
 *      * Trên khoảng đó nó song ánh = đơn điệu (x tăng thì y tăng, x giảm thì y giảm) + không bị lặp y??
 *              ! sign = (y_1-y_0)/(x_1 - x_0);
 *              ! left = 0; right = 1;
 *            * duyệt với i = 1,... với i <= số điểm (x,y)-2
 *                  ! signCur =(y_{i+1} - y_{i})/(x_{i+1} - x_{i})
 *                  * nếu signCurr!=sign thì
 *                      * Kết thúc một khoảng đơn điệu, lưu index của đầu mút bên phải: 
 *                          ! right = i;
 *                          ? Nếu sign != 0 (đủ điệu kiện là song ánh)
 *                              ! thêm (left, right) vào danh sách chứa các khoảng đơn điệu: ThisSpaces  
 *                      * Chuẩn bị sang khoảng mới
 *                          ! left = i;
 *                          ! sign = signCurr;
 *                      ? Nếu i == count-2 thì right = i+1 và thêm (left,right) vào khoảng             
 *                  ! i++;
 *      * Trên các khoảng đó có nghiệm ??
 *              * Duyệt các khoảng với i = 0,1,.. với i <= lastSpace = số khoảng(left,right)-1
 *                  ? Nghiệm nằm trong khoảng này không
 *                      ! nếu y không thuộc (y_left, y_right) hoặc (y_right, y_left) 
 *                              * Xóa khoảng này đi
 *                              * cập nhật lastSpace = lastSpace - 1;
 *                      ! nếu thuộc thì đây sẽ là một song ánh và có nghiệm 
 *                  ! i++;
 *      * Thu được các khoảng là song ánh và chứa nghiệm
 * 
 *  TODO-2: Từ các khoảng đc cho là song ánh => rút ra lượng mốc phù hợp để xây dựng nội suy
 *      * Với khoảng song ánh (left, right)
 *      * count = right - left + 1  = số mốc trong khoảng song ánh này
 *      * indexMinDist = left;
 *      * minDist = |_y - y_left|
 *      ? Nếu count <= countNumPoint 
 *            * thì lấy từ (left, right) để nội suy ngược
 *      ! Nếu không
 *          ! Điểm gần _y nhất tính từ left  
 *                 * Duyệt i = left+1, left+2,...,right
 *                      ! tempMin =  |_y - y_i|
 *                      ? tempMin < minDist: minDist = tempMin; indexMinDist = i; 
 *          ? Nếu (right - indexMinDist + 1 <= countNumPoint) 
 *                 * thì chọn (indexMinDist, right) để nội suy ngược (Mốc bất kì hoặc Newton Tiến)
 *          ? Nếu (indexMinDist - left + 1 <= countNumPoint) 
 *                 * thì chọn (left, indexMinDist) để nội suy ngược (Mốc bất kì hoặc Newton Lùi)
 *          ? Nếu không 
 *                 * thì chọn (indexMinDist, indexMinDist + counNumPoint - 1) để nội suy ngược (Mốc bất kì hoặc Newton Tiến)
 *     TODO: làm tương tự với các khoảng song ánh còn lại
 *          * Thu được danh sách các khoảng song ánh với lượng mốc <= countNumPoint để nội suy ngược
 *              
 * 
 *  TODO-3: Nội suy ngược trên khoảng song ánh đó
 *      * INPUT: (left, right)
 *      ? Nếu các mốc trên đây là cách đều 
 *          * Nội suy ngược với Newton tiến trên đây => return x;
 *      ? Nếu khongo
 *          * Nội suy ngược với Newton mốc bất kỳ hoặc Lagrange => return x;
 * 
 *  TODO: Nội suy ngược với mốc cách đều: (lặp Newton)
 *      * poly = Newtontien((left,right))
 *      * eta = poly - (y_0+t\Delta y_0)
 *      * t_0 = (y* - y_0)
 *      * t = t_0
 *      * x_curr = 0;
 *      * x_bef = 0;
 *      * eps = 99;
 *      ? Lặp t=\phi(t)
 *          * Khi eps >= epsilon
 *                  ! x_bef = x_0+h*t;
 *                  ! t = (t_0 - eta(t))/Delta y_0;
 *                  ! x_cur = x_0+h*t;
 *                  ! eps = |x_curr-x_bef|/x_curr
 *          * return x_curr;
 *          
 *  TODO: Nội suy ngược với mốc bất kỳ (Lagrange hoặc Newton)
**/