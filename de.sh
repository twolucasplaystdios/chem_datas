deactivate
{
    # 退出虛擬環境
    deactivate
} || {
    echo "退出虛擬環境時發生錯誤！錯誤代碼為： $?"
}

clear

echo "設置完成"
