#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <vector>
#include <string>
#include <thread>
#include <future>
#include <queue>
#include <atomic>
#include <functional>
#include <algorithm>

namespace torali {

  struct ThreadPool {
    bool stop;
    std::atomic<size_t> active;
    std::size_t maxQueueSize;
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    std::condition_variable condition_finished;
    std::condition_variable condition_producer;

    ThreadPool(size_t numThreads, std::size_t maxQ = 0) : stop(false), active(0), maxQueueSize(maxQ ? maxQ : numThreads * 4) {
      for (size_t i = 0; i < numThreads; ++i) {
	workers.emplace_back([this]() {
	  for (;;) {
	    std::function<void()> task;
	    {
	      std::unique_lock<std::mutex> lock(queue_mutex);
	      condition.wait(lock, [this]() { return stop || !tasks.empty(); });
	      if (stop && tasks.empty()) return;
	      task = std::move(tasks.front());
	      tasks.pop();
	      ++active;
	      condition_producer.notify_one();
	    }
	    task();
	    {
	      std::unique_lock<std::mutex> lock(queue_mutex);
	      --active;
	      condition_finished.notify_all();
	    }
	  }
	});
      }
    }

    template<class F>
    std::future<void> enqueue(F&& f) {
      auto task = std::make_shared<std::packaged_task<void()>>(std::forward<F>(f));
      std::future<void> res = task->get_future();
      {
	std::unique_lock<std::mutex> lock(queue_mutex);
	condition_producer.wait(lock, [this]() { return stop || tasks.size() < maxQueueSize; });
	tasks.emplace([task]() { (*task)(); });
      }
      condition.notify_one();
      return res;
    }

    void waitAll() {
      std::unique_lock<std::mutex> lock(queue_mutex);
      condition_finished.wait(lock, [this]() { return tasks.empty() && active == 0; });
    }

    ~ThreadPool() {
      {
	std::unique_lock<std::mutex> lock(queue_mutex);
	stop = true;
      }
      condition.notify_all();
      for (std::thread& worker : workers) worker.join();
    }
  };

}

#endif
