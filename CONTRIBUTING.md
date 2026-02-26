# Contributing Guide

感谢你为本项目做出贡献！本项目采用 **Fork + Pull Request** 协作模式。

组织成员为 Reporter 权限，不能直接 push 主仓库，必须通过 PR 合并。

---

## 🧭 开发流程

### 1. Fork 仓库
点击页面右上角 Fork 到个人账号。

---

### 2. 克隆 Fork

```bash
git clone https://github.com/<your_name>/<repo>.git
cd <repo>
```

---

### 3. 添加上游仓库（必须）

```bash
git remote add upstream https://github.com/<org>/<repo>.git
```

验证：
```bash
git remote -v
```

---

### 4. 同步上游最新变更

```bash
git fetch upstream
git checkout dev
git merge upstream/dev
git push origin dev
```

---

### 5. 创建功能分支（禁止在 main 和 dev 开发）

```bash
git checkout -b feature/<name>
```

命名规范：
```bash
feature/<功能名>
fix/<问题名>
refactor/<模块名>
docs/<文档名>
```

---

### 6. 开发并提交

```bash
git add .
git commit -m "feat: add xxx"
```

---

### 7. 推送到 Fork

```bash
git push origin feature/<name>
```

---

### 8. 创建 Pull Request

提交 PR 时必须包含：

- 功能说明
- 修改内容
- 测试方法
- 影响范围

---

### 9. 处理 Review

根据 reviewer 意见修改：

```bash
git add .
git commit -m "fix: address review comments"
git push origin feature/<name>
```

---

## 🔁 同步主仓库（建议每天执行）

```bash
git fetch upstream
git checkout dev
git merge upstream/dev
git push origin dev
```

---

## ⚠️ 重要规则

- 不允许直接 push 主仓库
- 所有代码必须通过 PR
- 每个 PR 只做一件事
- PR 尽量控制在 400 行以内
- 提交必须符合 commit 规范

---

## 🧪 冲突解决

```bash
git fetch upstream
git checkout feature/<name>
git merge upstream/dev
```

手动解决冲突后：
```bash
git add .
git commit -m "merge upstream main"
git push origin feature/<name>
```
---

感谢你的贡献 🚀