const deck = document.getElementById('deck');
const slides = Array.from(document.querySelectorAll('.slide'));
const outline = document.getElementById('outline');
const outlineToggle = document.getElementById('outlineToggle');
const outlineList = document.getElementById('outlineList');
const prevBtn = document.getElementById('prevBtn');
const nextBtn = document.getElementById('nextBtn');
const pageIndicator = document.getElementById('pageIndicator');
const progressBar = document.getElementById('progressBar');

let index = 0;
let wheelLock = false;
let touchStartY = null;

function renderOutline() {
  outlineList.innerHTML = '';
  slides.forEach((slide, i) => {
    const li = document.createElement('li');
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.textContent = `${String(i + 1).padStart(2, '0')} · ${slide.dataset.title || '未命名页'}`;
    btn.addEventListener('click', () => {
      go(i);
      outline.classList.remove('open');
    });
    if (i === index) btn.classList.add('active');
    li.appendChild(btn);
    outlineList.appendChild(li);
  });
}

function go(next) {
  const target = Math.max(0, Math.min(slides.length - 1, next));
  if (target === index) return;
  index = target;
  deck.style.transform = `translateY(-${index * 100}vh)`;
  pageIndicator.textContent = `${index + 1} / ${slides.length}`;
  progressBar.style.width = `${((index + 1) / slides.length) * 100}%`;
  renderOutline();
}

function initState() {
  deck.style.transform = 'translateY(0)';
  pageIndicator.textContent = `1 / ${slides.length}`;
  progressBar.style.width = `${(1 / slides.length) * 100}%`;
  renderOutline();
}

outlineToggle.addEventListener('click', () => outline.classList.toggle('open'));
prevBtn.addEventListener('click', () => go(index - 1));
nextBtn.addEventListener('click', () => go(index + 1));

document.addEventListener('keydown', (e) => {
  if (['ArrowRight', 'ArrowDown', 'PageDown', ' '].includes(e.key)) {
    e.preventDefault();
    go(index + 1);
  }
  if (['ArrowLeft', 'ArrowUp', 'PageUp'].includes(e.key)) {
    e.preventDefault();
    go(index - 1);
  }
  if (e.key.toLowerCase() === 'o') outline.classList.toggle('open');
  if (e.key === 'Escape') outline.classList.remove('open');
});

document.addEventListener('wheel', (e) => {
  if (wheelLock || Math.abs(e.deltaY) < 28) return;
  wheelLock = true;
  go(index + (e.deltaY > 0 ? 1 : -1));
  setTimeout(() => { wheelLock = false; }, 420);
}, { passive: true });

document.addEventListener('touchstart', (e) => {
  touchStartY = e.changedTouches[0].clientY;
}, { passive: true });

document.addEventListener('touchend', (e) => {
  if (touchStartY === null) return;
  const endY = e.changedTouches[0].clientY;
  const delta = touchStartY - endY;
  if (Math.abs(delta) > 40) go(index + (delta > 0 ? 1 : -1));
  touchStartY = null;
}, { passive: true });

Array.from(document.querySelectorAll('img')).forEach((img) => {
  img.addEventListener('error', () => {
    img.classList.add('is-broken');
    if (!img.closest('figure')?.querySelector('.img-error')) {
      const p = document.createElement('p');
      p.className = 'img-error';
      p.textContent = '图片加载失败';
      p.style.cssText = 'position:absolute;inset:auto 8px 8px 8px;margin:0;padding:6px 8px;background:rgba(180,40,40,.75);color:#fff;border-radius:6px;font-size:12px;z-index:2';
      img.closest('figure')?.appendChild(p);
    }
  });
});

initState();
