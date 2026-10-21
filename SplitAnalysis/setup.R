# setup.R


options(repos = c(CRAN = "https://cloud.r-project.org"))


if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("错误: 'remotes' 包未找到。请确保 environment.yaml 已正确安装。")
}


cat("正在检查编译器环境...\n")

system("which gcc") 


pkg_path <- "./EccDNAFeature"

if (file.exists(pkg_path)) {
  message("正在安装本地包: ", pkg_path)
  message("依赖关系已通过 Conda 托管，R 将跳过依赖安装...")
  
  remotes::install_local(
    pkg_path, 
    dependencies = FALSE,  
    upgrade = "never",    
    build_vignettes = FALSE 
  )
  
  message("\n安装完成！请尝试运行 library(EccDNAFeature) 验证。")
  
} else {
  stop(paste("错误: 找不到目录", pkg_path, "。请确保您在包含该文件夹的目录下运行此脚本。"))
}