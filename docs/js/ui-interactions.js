document.addEventListener('DOMContentLoaded', () => {
    // 1. Scroll Reveal Animation
    const observerOptions = {
        root: null,
        rootMargin: '0px',
        threshold: 0.1
    };

    const observer = new IntersectionObserver((entries, observer) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.classList.add('visible');
                observer.unobserve(entry.target); // Only animate once
            }
        });
    }, observerOptions);

    const scrollElements = document.querySelectorAll('.animate-on-scroll');
    scrollElements.forEach(el => observer.observe(el));

    // 2. Parallax Effect for Hero Section (Disabled)
    /*
    const heroSection = document.querySelector('.hero-section');
    const heroTitle = document.querySelector('.hero-title');
    const heroSubtitle = document.querySelector('.hero-subtitle');
    
    if (heroSection && heroTitle && heroSubtitle) {
        heroSection.addEventListener('mousemove', (e) => {
            const x = (window.innerWidth - e.pageX * 2) / 100;
            const y = (window.innerHeight - e.pageY * 2) / 100;

            heroTitle.style.transform = `translate(${x * 0.5}px, ${y * 0.5}px)`;
            heroSubtitle.style.transform = `translate(${x * 0.2}px, ${y * 0.2}px)`;
        });
    }
    */

    // 3. 3D Tilt Effect for Cards (Event Delegation)
    document.addEventListener('mousemove', (e) => {
        // Excluded .comm-card and .citation-card-new as per request
        const card = e.target.closest('.feature-item, .doc-card, .pub-item');
        if (!card) return;

        const rect = card.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;
        
        const centerX = rect.width / 2;
        const centerY = rect.height / 2;
        
        const rotateX = ((y - centerY) / centerY) * -5; // Max 5deg rotation
        const rotateY = ((x - centerX) / centerX) * 5;

        card.style.transform = `perspective(1000px) rotateX(${rotateX}deg) rotateY(${rotateY}deg) scale3d(1.02, 1.02, 1.02)`;
    });

    document.addEventListener('mouseout', (e) => {
        // Excluded .comm-card and .citation-card-new as per request
        const card = e.target.closest('.feature-item, .doc-card, .pub-item');
        if (!card) return;
        
        // Reset transform when mouse leaves the element
        // Note: mouseout bubbles, so we need to check if we really left the card
        if (!card.contains(e.relatedTarget)) {
            card.style.transform = 'perspective(1000px) rotateX(0) rotateY(0) scale3d(1, 1, 1)';
        }
    });
});
